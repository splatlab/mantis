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

using LRUCacheMap =  LRU::Cache<uint64_t, std::vector<uint64_t>>;
using SpinLockT = std::mutex;

typedef sdsl::bit_vector BitVector;
typedef sdsl::rrr_vector<63> BitVectorRRR;
typedef uint32_t colorIdType;


struct Edge {
    colorIdType n1;
    colorIdType n2;

    Edge() : n1{static_cast<colorIdType>(-1)}, n2{static_cast<colorIdType>(-1)} {}
    Edge(colorIdType inN1, colorIdType inN2) : n1(inN1), n2(inN2) {}

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

class MST {
public:
    MST(std::string prefix, std::shared_ptr<spdlog::logger> logger, uint32_t numThreads);

    void buildMST();

private:
    bool buildEdgeSets();

    void findNeighborEdges(CQF<KeyObject> &cqf, KeyObject &keyobj, std::vector<Edge> &edgeList);

    bool calculateWeights();

    bool encodeColorClassUsingMST();

    DisjointSets kruskalMSF();

    std::set<workItem> neighbors(CQF<KeyObject> &cqf, workItem n);

    bool exists(CQF<KeyObject> &cqf, dna::canonical_kmer e, uint64_t &eqid);

    uint64_t hammingDist(uint64_t eqid1, uint64_t eqid2,
                         uint64_t &srcId, std::vector<uint64_t> &srcEq);

    std::vector<uint32_t> getDeltaList(uint64_t eqid1, uint64_t eqid2);

    void buildColor(std::vector<uint64_t> &eq, uint64_t eqid, BitVectorRRR *bv);

    inline uint64_t getBucketId(uint64_t c1, uint64_t c2);

    void buildPairedColorIdEdgesInParallel(uint32_t threadId, CQF<KeyObject> &cqf,
                                           std::vector<spp::sparse_hash_set<Edge, edge_hash>> &edgesetList,
                                           sdsl::bit_vector &nodes, uint64_t &maxId, uint64_t &numOfKmers);

    void calcHammingDistInParallel(uint32_t i, std::vector<Edge> &edgeList);

    void calcDeltasInParallel(uint32_t threadID, uint64_t cbvID1, uint64_t cbvID2,
            sdsl::int_vector<> &parentbv, sdsl::int_vector<> &deltabv,
            sdsl::bit_vector::select_1_type &sbbv );

    std::string prefix;
    uint32_t numSamples = 0;
    uint64_t k;
    uint64_t num_of_ccBuffers;
    uint64_t num_edges = 0;
    uint64_t num_colorClasses = 0;
    uint64_t mstTotalWeight = 0;
    colorIdType zero = static_cast<colorIdType>(UINT64_MAX);
    BitVectorRRR *bvp1, *bvp2;
    LRUCacheMap lru_cache;
    uint64_t gcntr = 0;
    std::vector<std::string> eqclass_files;
    std::vector<std::vector<Edge>> edgeBucketList;
    std::vector<std::vector<Edge>> weightBuckets;
    std::vector<std::vector<std::pair<colorIdType, uint32_t> >> mst;
    spdlog::logger *logger{nullptr};
    uint32_t nThreads = 1;
    SpinLockT colorMutex;

};

#endif //MANTIS_MST_H
