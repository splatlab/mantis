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

using FilterType = CQF<KeyObject>;

typedef uint32_t colorIdType;
typedef uint32_t weightType;
typedef unsigned __int128 uint128_t;

struct Cost {
    uint64_t numSteps{0};
    uint64_t numQueries{0};
};

struct CompactEdge {

    static uint128_t getEdge(uint64_t hi, uint64_t lo, uint64_t nodeSize) {
        return (uint128_t)(lo) + (((uint128_t)hi) << nodeSize);
    }

    static std::pair<uint64_t, uint64_t> getNodes(uint128_t edge, uint64_t nodeSize) {
        uint64_t mask = ((1ULL << nodeSize) - 1);
        uint64_t hi = (edge >> nodeSize) & mask;
        uint64_t lo = edge & mask;
        return std::make_pair(hi, lo);
    }
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
/*
// note: @fatal: careful! The hash highly depends on the length of the edge ID (uint32)
struct edge_hash {
    uint64_t operator()(const Edge &e) const {
        return MurmurHash64A(&e, sizeof(Edge), 2038074743);
        *//*uint64_t res = e.n1;
        return (res << 32) | (uint64_t)e.n2;*//*
    }
};*/

struct workItem {
    dna::canonical_kmer node;
    colorIdType colorId;

    workItem(dna::canonical_kmer n, colorIdType c) : node(n), colorId(c) {}

    // Required to be able to use it as a key in set
    bool operator<(const workItem &item2) const {
        return (*this).node < item2.node;
    }
};

/*

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
*/

// To represent Disjoint Sets
struct DisjointSets {
    sdsl::int_vector<> els;
//    std::vector<DisjointSetNode> els; // the size of total number of colors
//    uint64_t n;

    // Constructor.
    explicit DisjointSets(uint64_t n) {
        // Allocate memory
//        this->n = n;
        els = sdsl::int_vector<>(n, 1, ceil(log2(n))+1); // 1 bit which is set if the node IS its own parent
//        els.resize(n);
        // Initially, all vertices are in
        // different sets and have rank 0.
//        for (uint64_t i = 0; i < n; i++) {
//            //every element is parent of itself
//            els[i].setParent(static_cast<colorIdType>(i));
//        }
    }

    bool selfParent(uint64_t idx) {
        if (idx >= els.size()) {
            std::cerr << "ERROR in selfParent => idx > vector size: " << idx << " " << els.size() << "\n";
            std::exit(3);
        }
        return els[idx] & 1ULL;
    }

    void setParent(uint64_t idx, uint64_t parentIdx) {
        bool ownParent = idx == parentIdx;
        if (idx >= els.size() or parentIdx >= els.size()) {
            std::cerr << "ERROR in setParent => idx > vector size: "
            << idx << " " << parentIdx << " " << els.size() << "\n";
        }
        parentIdx = getParent(parentIdx);
        els[idx] = (parentIdx << 1ULL) | static_cast<uint64_t>(ownParent);
    }

    uint64_t getParent(uint64_t idx) {
        if (idx >= els.size()) {
            std::cerr << "ERROR in getParent => idx > vector size: " << idx << " " << els.size() << "\n";
            std::exit(3);
        }
        if (selfParent(idx))
            return idx;
        return els[idx] >> 1ULL;
    }

    uint64_t getRank(uint64_t idx) {
        if (idx >= els.size()) {
            std::cerr << "ERROR in getRank => idx > vector size: " << idx << " " << els.size() << "\n";
            std::exit(3);
        }
        uint64_t par = find(idx);
        return els[par] >> 1ULL;
    }

    void incrementRank(uint64_t idx) {
        if (idx >= els.size()) {
            std::cerr << "ERROR in incrementRank => idx > vector size: " << idx << " " << els.size() << "\n";
            std::exit(3);
        }
        auto parIdx = find(idx);
        uint64_t rank = els[parIdx] >> 1ULL;
        ++rank;
        els[parIdx] = (rank << 1ULL) | 1ULL; // selfParent bit is set
    }
    // Find the parent of a node 'u'
    // Path Compression
    uint64_t find(uint64_t u) {
        /* Make the parent of the nodes in the path
           from u--> parent[u] point to parent[u] */
        if (not selfParent(u)) {
            setParent(u, find(getParent(u)));
        }
        return getParent(u);
    }

    // Union by rank
    void merge(uint64_t x, uint64_t y, uint32_t edgeW) {
        auto parent = find(x), child = find(y);

        /* Make tree with smaller height
           a subtree of the other tree  */
        if (getRank(child) > getRank(parent)) {
            std::swap(child, parent);
        }
        if (getRank(child) == getRank(parent)) {
            incrementRank(parent);
        }
        setParent(child, parent);
    }
};

/**
 * Adjacency List
 * Weight is stored within the smallerSrc vector
 */
struct AdjList {
    // weight is in the smallerSrc
    sdsl::int_vector<> smallerSrc;
    sdsl::bit_vector smallerSrcBoundary;
    sdsl::bit_vector::select_1_type smallerSrcBoundarySel;
    sdsl::int_vector<> greaterSrc;
    sdsl::bit_vector greaterSrcBoundary;
    sdsl::bit_vector::select_1_type greaterSrcBoundarySel;

    uint64_t weightBits;
    uint64_t weightMask;

    AdjList(uint64_t numColorClasses, uint64_t numSamples) {
        weightBits = static_cast<uint64_t>(ceil(log2(numSamples)));
        weightMask = (1ULL << weightBits) - 1;
        smallerSrc = sdsl::int_vector<>(numColorClasses, 0, ceil(log2(numColorClasses))+weightBits);
        smallerSrcBoundary = sdsl::bit_vector(numColorClasses);
        greaterSrc = sdsl::int_vector<>(numColorClasses, 0, ceil(log2(numColorClasses))+1);
        greaterSrcBoundary = sdsl::bit_vector(numColorClasses);
    }

    bool loadFromFile(const char* file) {
        std::ifstream input(file);
        uint64_t src,dest, weight, smallerIdx{0}, greaterIdx{0}, prevSrc{0}, smallerNeighborCnt{0}, greaterNeighborCnt{0};
        input >> src >> dest >> weight;
        while (input.good()) {
            if (smallerIdx == smallerSrc.size() or greaterIdx == greaterSrc.size()) {
                std::cerr << "ERROR! Found more edges than expected for the MST stored in the temp file " << file << "\n";
                std::exit(3);
            }
            if (src < dest) {
                if (prevSrc != src) {
                    smallerIdx += (smallerNeighborCnt == 0); // leave the slot zero if no greater neighbor nodes
                    smallerSrcBoundary[smallerIdx] = 1;
                    smallerNeighborCnt = 0;
                }
                if (smallerIdx == smallerSrc.size() ) {
                    std::cerr << "ERROR! src<dest: Found more edges than expected for the MST stored in the temp file " << file << "\n";
                    std::exit(3);
                }
                smallerSrc[smallerIdx] = (dest << weightBits) & weight;
                ++smallerNeighborCnt;
                ++smallerIdx;
            } else {
                if (prevSrc != src) {
                    greaterIdx += (greaterNeighborCnt == 0); // leave the slot zero if no smaller neighbor nodes
                    greaterSrcBoundary[greaterIdx] = 1;
                    greaterNeighborCnt = 0;
                }
                if (greaterIdx == greaterSrc.size() ) {
                    std::cerr << "ERROR! src>dest: Found more edges than expected for the MST stored in the temp file " << file << "\n";
                    std::exit(3);
                }
                greaterSrc[greaterIdx] = (dest << 1ULL) & 1ULL;
                ++greaterNeighborCnt;
                ++greaterIdx;
            }
            input >> src >> dest >> weight;
        }
        input.close();
        smallerSrcBoundarySel = sdsl::bit_vector::select_1_type(&smallerSrcBoundary);
        greaterSrcBoundarySel = sdsl::bit_vector::select_1_type(&greaterSrcBoundary);
        return greaterIdx == greaterSrc.size() and smallerIdx == smallerSrc.size();
    }

    std::vector<colorIdType> neighbors(uint64_t idx) {
        std::vector<colorIdType> res;
        res.reserve(1000);
        auto start = smallerSrcBoundarySel(idx);
        bool found = false;
        uint64_t wrd{0};
        do {
            uint64_t wrdLen = std::min(static_cast<uint64_t >(64), smallerSrcBoundary.size()-start);
            wrd = smallerSrcBoundary.get_int(start, wrdLen);
            for (uint64_t j = 0; j < wrdLen; j++) {
                auto nei = smallerSrc[start + j];
                if (nei & weightMask)
                    res.push_back(nei >> weightBits);
                if ((wrd >> j) & 0x01) {
                    found = true;
                    break;
                }
            }
            start += wrdLen;
        } while (!found);
    //// greater
        start = greaterSrcBoundarySel(idx);
        found = false;
        do {
            uint64_t wrdLen = std::min(static_cast<uint64_t >(64), greaterSrcBoundary.size()-start);
            wrd = greaterSrcBoundary.get_int(start, wrdLen);
            for (uint64_t j = 0; j < wrdLen; j++) {
                auto nei = greaterSrc[start + j];
                if (nei & 1ULL)
                    res.push_back(nei >> 1ULL);
                if ((wrd >> j) & 0x01) {
                    found = true;
                    break;
                }
            }
            start += wrdLen;
        } while (!found);

        return res;
    }

    uint64_t getWeight(uint64_t n1, uint64_t n2) {
        if (n1 == n2)
            return 0;
        if (n1 > n2) {
            std::swap(n1, n2);
        }
        auto start = smallerSrcBoundarySel(n1);
        bool found = false;
        uint64_t wrd{0};
        do {
            uint64_t wrdLen = std::min(static_cast<uint64_t >(64), smallerSrcBoundary.size()-start);
            wrd = smallerSrcBoundary.get_int(start, wrdLen);
            for (uint64_t j = 0; j < wrdLen; j++) {
                auto nei = smallerSrc[start + j];
                if (nei & weightMask) {
                    auto neighbor = (nei >> weightBits);
                    if (neighbor == n2) {
                        return (nei & weightMask);
                    }
                }
                if ((wrd >> j) & 0x01) {
                    found = true;
                    break;
                }
            }
            start += wrdLen;
        } while (!found);
        return invalid;
    }

};

class MSTMerger {
public:
    MSTMerger(std::string prefix,
            spdlog::logger *logger,
            uint32_t numThreads,
            std::string prefix1,
            std::string prefix2);

    void mergeMSTs();

private:

    std::pair<uint64_t, uint64_t> buildMultiEdgesFromCQFs();

    bool buildEdgeSets();

    bool calculateMSTBasedWeights();

    bool encodeColorClassUsingMST();

    void kruskalMSF();

    std::set<workItem> neighbors(CQF<KeyObject> &cqf, workItem n);

    bool exists(CQF<KeyObject> &cqf, dna::canonical_kmer e, uint64_t &eqid);

    uint64_t mstBasedHammingDist(uint64_t eqid1,
                                 uint64_t eqid2,
                                 MSTQuery *mst,
                                 LRUCacheMap &lru_cache,
                                 std::vector<uint64_t> &srcEq,
                                 QueryStats &queryStats,
                                 std::unordered_map<uint64_t, std::vector<uint64_t>> &fixed_cache);

    void buildPairedColorIdEdgesInParallel(uint32_t threadId, CQF<KeyObject> &cqf,
                                           uint64_t &cnt, uint64_t &maxId, uint64_t &numOfKmers);


    void findNeighborEdges(CQF<KeyObject> &cqf, KeyObject &keyobj, std::vector<Edge> &edgeList);

    void calcMSTHammingDistInParallel(uint32_t i,
                                  std::vector<std::pair<colorIdType, weightType>> &edgeList,
                                  std::vector<uint32_t> &srcStarts,
                                  MSTQuery *mst,
                                  std::vector<LRUCacheMap> &lru_cache,
                                  std::vector<QueryStats> &queryStats,
                                  std::unordered_map<uint64_t, std::vector<uint64_t>> &fixed_cache,
                                  uint32_t numSamples);

    void calcDeltasInParallel(uint32_t threadID, uint64_t s, uint64_t e, uint64_t deltaOffset,
                              sdsl::int_vector<> &parentbv, sdsl::int_vector<> &deltabv,
                              sdsl::bit_vector &bbv,
                              bool isMSTBased);

    void buildMSTBasedColor(uint64_t eqid1, MSTQuery *mst1,
                            LRUCacheMap &lru_cache1, std::vector<uint64_t> &eq1,
                            QueryStats &queryStats,
                            std::unordered_map<uint64_t, std::vector<uint64_t>> &fixed_cache);

    std::vector<uint32_t> getMSTBasedDeltaList(uint64_t eqid1, uint64_t eqid2,
                                               MSTQuery * mstPtr,
                                               std::unordered_map<uint64_t, std::vector<uint64_t>>& fixed_cache,
                                               LRUCacheMap &lru_cache,
                                               QueryStats &queryStats);

    void planCaching(MSTQuery *mst,
                     std::vector<std::pair<colorIdType, weightType>> &edges,
                     std::vector<uint32_t> &srcStartIdx,
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
    uint64_t numCCPerBuffer;
    uint64_t num_edges = 0;
    uint64_t num_colorClasses = 0;
    uint64_t mstTotalWeight = 0;
    colorIdType zero = static_cast<colorIdType>(UINT64_MAX);
    std::vector<LRUCacheMap> lru_cache1;//10000);
    std::vector<LRUCacheMap> lru_cache2;//10000);
    std::unordered_map<uint64_t, std::vector<uint64_t>> fixed_cache1;
    std::unordered_map<uint64_t, std::vector<uint64_t>> fixed_cache2;
    std::vector<QueryStats> queryStats1;
    std::vector<QueryStats> queryStats2;
    std::vector<std::pair<colorIdType, colorIdType >> colorPairs;
    std::string prefix1;
    std::string prefix2;
    std::unique_ptr<MSTQuery> mst1;
    std::unique_ptr<MSTQuery> mst2;
    std::unique_ptr<std::vector<Edge>> edges;
    std::vector<std::unique_ptr<std::vector<Edge>>> weightBuckets;
    std::unique_ptr<sdsl::int_vector<>> mst;
    std::unique_ptr<sdsl::bit_vector> mstBbv;
//    std::unique_ptr<std::vector<std::vector<std::pair<colorIdType, uint32_t> >>> mst;
    spdlog::logger *logger{nullptr};
    uint32_t nThreads = 1;
    SpinLockT colorMutex;

    uint64_t numBlocks;

};

#endif //MANTIS_MST_H
