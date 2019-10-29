//
// Created by Fatemeh Almodaresi on 2018-10-13.
//

#ifndef MANTIS_MST_H
#define MANTIS_MST_H

#include <set>
#include <vector>
#include <queue>
#include <cstdint>
#include <mutex>
#include <thread>

// sparsepp should be included before gqf_cpp! ow, we'll get a conflict in MAGIC_NUMBER
#include "sparsepp/spp.h"

#include "mantisconfig.hpp"
#include "gqf_cpp.h"
#include "spdlog/spdlog.h"

#include "canonicalKmer.h"
#include "sdsl/bit_vectors.hpp"
#include "gqf/hashutil.h"

#include "lru/lru.hpp"
#include "mstQuery.h"

using LRUCacheMap =  LRU::Cache<uint64_t, std::vector<uint64_t>>;

//using LRUCacheMap = HPHP::ConcurrentScalableCache<uint64_t , std::vector<uint64_t >>;

using SpinLockT = std::mutex;

typedef sdsl::bit_vector BitVector;
typedef sdsl::rrr_vector<63> BitVectorRRR;
typedef uint32_t colorIdType;
typedef uint16_t weightType;

struct Cost {
    uint64_t numSteps{0};
    uint64_t numQueries{0};
};

// undirected edge
struct Edge {
    colorIdType n1;
    colorIdType n2;

    Edge() : n1{static_cast<colorIdType>(-1)}, n2{static_cast<colorIdType>(-1)} {}

    Edge(colorIdType inN1, colorIdType inN2) : n1(inN1), n2(inN2) {
        if (n1 > n2) {
            std::swap(n1, n2);
        }
    }

    bool operator==(const Edge &e) const {
        return n1 == e.n1 && n2 == e.n2;
    }
};

// note: @fatal: careful! The hash highly depends on the length of the edge ID (uint32)
struct edge_hash {
    uint64_t operator()(const Edge &e) const {
        return MurmurHash64A(&e, sizeof(Edge), 2038074743);
        /*uint64_t res = e.n1;
        return (res << 32) | (uint64_t)e.n2;*/
    }
};

struct workItem {
    dna::canonical_kmer node;
    colorIdType colorId;

    workItem(dna::canonical_kmer n, colorIdType c) : node(n), colorId(c) {}

    // Required to be able to use it as a key in set
    bool operator<(const workItem &item2) const {
        return (*this).node < item2.node;
    }
};


struct DisjointSetNode {
    colorIdType parent{0};
    uint64_t rnk{0}, w{0}, edges{0};

    void setParent(colorIdType p) { parent = p; }

    void mergeWith(DisjointSetNode &n, uint32_t edgeW, colorIdType id) {
        n.setParent(parent);
        w += (n.w + static_cast<uint64_t>(edgeW));
        edges += (n.edges + 1);
        n.edges = 0;
        n.w = 0;
        if (rnk == n.rnk) {
            rnk++;
        }
    }
};

// To represent Disjoint Sets
struct DisjointSets {
    std::vector<DisjointSetNode> els;
    uint64_t n;

    // Constructor.
    explicit DisjointSets(uint64_t n) {
        // Allocate memory
        this->n = n;
        els.resize(n);
        // Initially, all vertices are in
        // different sets and have rank 0.
        for (uint64_t i = 0; i < n; i++) {
            //every element is parent of itself
            els[i].setParent(static_cast<colorIdType>(i));
        }
    }

    // Find the parent of a node 'u'
    // Path Compression
    uint32_t find(colorIdType u) {
        /* Make the parent of the nodes in the path
           from u--> parent[u] point to parent[u] */
        if (u != els[u].parent)
            els[u].parent = find(els[u].parent);
        return els[u].parent;
    }

    // Union by rank
    void merge(colorIdType x, colorIdType y, uint32_t edgeW) {
        x = find(x), y = find(y);

        /* Make tree with smaller height
           a subtree of the other tree  */
        if (els[x].rnk > els[y].rnk) {
            els[x].mergeWith(els[y], edgeW, x);

        } else {// If rnk[x] <= rnk[y]
            els[y].mergeWith(els[x], edgeW, y);
        }
    }
};

class MSTMerger {
public:
    MSTMerger(CQF<KeyObject> *cqf, std::string prefix, spdlog::logger *logger, uint32_t numThreads, std::string prefix1,
              std::string prefix2, uint64_t numColorBuffers = 1);

    void mergeMSTs();

private:
    bool buildEdgeSets();

    void findNeighborEdges(CQF<KeyObject> &cqf, KeyObject &keyobj, std::vector<Edge> &edgeList);

    bool calculateMSTBasedWeights();

    bool encodeColorClassUsingMST();

    DisjointSets kruskalMSF();

    std::set<workItem> neighbors(CQF<KeyObject> &cqf, workItem n);

    bool exists(CQF<KeyObject> &cqf, dna::canonical_kmer e, uint64_t &eqid);

    uint64_t mstBasedHammingDist(uint64_t eqid1,
                                 uint64_t eqid2,
                                 MSTQuery *mst,
                                 LRUCacheMap &lru_cache,
                                 std::vector<uint64_t> &srcEq,
                                 QueryStats &queryStats,
                                 std::unordered_map<uint64_t, std::vector<uint64_t>> &fixed_cache);

    inline uint64_t getBucketId(uint64_t c1, uint64_t c2);

    void writePotentialColorIdEdgesInParallel(uint32_t threadId,
                                              CQF<KeyObject> &cqf,
                                              std::queue<uint64_t> &blockIds,
                                              std::vector<std::ofstream> &blockFiles,
                                              std::vector<uint64_t> &blockCnt);

    void buildPairedColorIdEdgesInParallel(uint32_t threadId,
                                           CQF<KeyObject> &cqf,
                                           std::queue<uint64_t> &blockIds,
                                           std::vector<std::ifstream> &blockFiles,
                                           std::vector<uint64_t> &blockCnt,
                                           uint64_t &maxId);

    void calcMSTHammingDistInParallel(uint32_t i,
                                      std::vector<std::pair<colorIdType, weightType>> &edgeList,
                                      std::vector<uint32_t> &srcStarts,
                                      MSTQuery *mst,
                                      std::vector<LRUCacheMap> &lru_cache,
                                      std::vector<QueryStats> &queryStats,
                                      std::unordered_map<uint64_t, std::vector<uint64_t>> &fixed_cache,
                                      uint32_t numSamples);

    void calcDeltasInParallel(uint32_t threadID, uint64_t cbvID1, uint64_t cbvID2,
                              sdsl::int_vector<> &parentbv, sdsl::int_vector<> &deltabv,
                              sdsl::bit_vector::select_1_type &sbbv,
                              bool isMSTBased);

    void buildMSTBasedColor(uint64_t eqid1, MSTQuery *mst1,
                            LRUCacheMap &lru_cache1, std::vector<uint64_t> &eq1,
                            QueryStats &queryStats,
                            std::unordered_map<uint64_t, std::vector<uint64_t>> &fixed_cache);

    std::vector<uint32_t> getMSTBasedDeltaList(uint64_t eqid1, uint64_t eqid2, LRUCacheMap &lru_cache,
                                               bool isFirst, QueryStats &queryStats);

    void planCaching(MSTQuery *mst,
                     std::vector<std::pair<colorIdType, weightType>> &edges,
                     std::vector<uint32_t> &outDegrees,
                     std::vector<colorIdType> &colorsInCache);

    void planRecursively(uint64_t nodeId,
                         std::vector<std::vector<colorIdType>> &children,
                         std::vector<Cost> &mstCost,
                         std::vector<colorIdType> &colorsInCache,
                         uint64_t &cntr);

    std::string prefix;
    uint32_t numSamples = 0;
    uint32_t numOfFirstMantisSamples = 0;
    uint32_t secondMantisSamples = 0;
    uint64_t k;
    uint64_t num_of_ccBuffers = 1;
    uint64_t num_edges = 0;
    uint64_t num_colorClasses = 0;
    uint64_t mstTotalWeight = 0;
    colorIdType zero = static_cast<colorIdType>(UINT64_MAX);
    BitVectorRRR *bvp1, *bvp2;
    LRUCacheMap lru_cache;//10000);
    std::vector<LRUCacheMap> lru_cache1;//10000);
    std::vector<LRUCacheMap> lru_cache2;//10000);
    uint64_t fixed_size = 200;
    std::unordered_map<uint64_t, std::vector<uint64_t>> fixed_cache1;
    std::unordered_map<uint64_t, std::vector<uint64_t>> fixed_cache2;
    std::vector<QueryStats> queryStats1;
    std::vector<QueryStats> queryStats2;
    uint64_t gcntr = 0;
    std::vector<std::string> eqclass_files;
    std::vector<std::pair<colorIdType, colorIdType >> colorPairs;
    std::string prefix1;
    std::string prefix2;
    MSTQuery *mst1;
    MSTQuery *mst2;
    CQF<KeyObject> *cqf;
    std::vector<std::vector<Edge>> edgeBucketList;
    std::vector<std::vector<Edge>> weightBuckets;
    std::vector<std::vector<std::pair<colorIdType, uint32_t> >> mst;
    spdlog::logger *logger{nullptr};
    uint32_t nThreads = 1;
    SpinLockT colorMutex;
    std::vector<SpinLockT*> blockMutex;

    std::unordered_set<uint64_t> test;

};

#endif //MANTIS_MST_H
