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
    sdsl::int_vector<> smallerSrcCnt;
    sdsl::int_vector<> greaterSrc;
    sdsl::int_vector<> greaterSrcCnt;

    uint64_t weightBits;
    uint64_t weightMask;
    uint64_t idx{0};

    AdjList(uint64_t numColorClasses, uint64_t numSamples) {
        weightBits = static_cast<uint64_t>(ceil(log2(numSamples)));
        weightMask = (1ULL << weightBits) - 1;
        smallerSrc = sdsl::int_vector<>(numColorClasses, 0, ceil(log2(numColorClasses))+weightBits);
        smallerSrcCnt = sdsl::int_vector<>(numColorClasses, 0, ceil(log2(numColorClasses)));
        greaterSrc = sdsl::int_vector<>(numColorClasses, 0, ceil(log2(numColorClasses)));
        greaterSrcCnt = sdsl::int_vector<>(numColorClasses, 0, ceil(log2(numColorClasses)));
    }

    /**
     * IMportant NOTE:: Assumes src < dest!!
     * @param src smaller vertex of the edge
     * @param dest greater vertex of the edge
     * @param weight of the edge
     */
    void addEdge(uint64_t src, uint64_t dest, uint64_t weight) {
        if (src > dest) {
            std::cerr << "WARNING!! Expect src < dest but " << src << " > " << dest << "\n";
            std::swap(src, dest);
        }
        smallerSrc[idx] = (((uint128_t)dest) << weightBits) | weight;
        smallerSrcCnt[src] = smallerSrcCnt[src]+1;
        greaterSrcCnt[dest] = greaterSrcCnt[dest]+1;
        ++idx;
    }

    void fillGreater() {
        std::cerr << "\n\nIDX\n\n" << idx << "\n";
        for (auto i = 1; i < smallerSrcCnt.size(); i++) {
            smallerSrcCnt[i] = smallerSrcCnt[i-1] + smallerSrcCnt[i];
            greaterSrcCnt[i] = greaterSrcCnt[i-1] + greaterSrcCnt[i];
        }
        std::cerr << "\n";
        uint64_t src = 0;
        uint64_t edgeCntr = 0;
        for (auto i = 0; i < smallerSrcCnt.size(); i++) {
            while (edgeCntr < smallerSrcCnt[i]) {
                uint64_t dest = (smallerSrc[edgeCntr] >> weightBits);
                greaterSrcCnt[dest] = greaterSrcCnt[dest] - 1;
                greaterSrc[greaterSrcCnt[dest]] = dest;
                edgeCntr++;
            }
        }
        for (auto i = smallerSrcCnt.size()-1; i > 0; i--) {
            smallerSrcCnt[i] = smallerSrcCnt[i-1];
        }
        smallerSrcCnt[0] = 0;
    }

    /**
     * dfs + bfs
     */
     void hybridTreeWalk(uint64_t root, sdsl::int_vector<> &parentbv, sdsl::bit_vector &visited, uint64_t level=0) {
        sdsl::bit_vector currentLevel(parentbv.size(), 0);
        uint64_t remaining = 0;
        parentbv[root] = root; // and it's its own parent (has no parent)
        uint64_t parIdx = root;
        dfs(parIdx, parentbv, visited, currentLevel, remaining, 0);
        while (remaining) {
            uint64_t idx = 0;
            while (remaining and idx < currentLevel.size()) {
                auto wrdLen = std::min(currentLevel.size()-idx, static_cast<uint64_t >(64));
                auto wrd = currentLevel.get_int(idx, wrdLen);
                if (sdsl::bits::cnt(wrd) != 0) {
                    for (auto i = 0; i < wrdLen and remaining; i++) {
                        parIdx = idx + i;
                        if ( ((wrd >> i) & 1ULL) and visited[parIdx] == 0) {
                            remaining--;
                            dfs(parIdx, parentbv, visited, currentLevel, remaining, 0);
                        }
                    }
                }
                auto newWrd = currentLevel.get_int(idx, wrdLen);
                currentLevel.set_int(idx, newWrd xor wrd, wrdLen); // this is so cool and smart :D
                idx += wrdLen;
            }
        }
     }
    /**
     * DFS walk instead of BFS
     */
    void dfs(uint64_t parIdx,
            sdsl::int_vector<> &parentbv,
            sdsl::bit_vector &visited,
            sdsl::bit_vector &activeLevel,
            uint64_t &remaining,
            uint64_t level=0) {
        if (level == 200) {
            activeLevel[parIdx] = 1;
            remaining++;
//            std::cerr << "remaining " << remaining << "\n";
            return;
        }
//        std::cerr << "l" << level++ << "  ===>   ";
//        std::cerr << "p" << parIdx << " ";
        if (parIdx >= visited.size()) {
            std::cerr << "1 happened\n";
            std::exit(3);
        }
        visited[parIdx] = 1;

        if (parIdx >= smallerSrcCnt.size()) {
            std::cerr << "2 happened\n";
            std::exit(3);
        }
        uint64_t startIdx = smallerSrcCnt[parIdx];
        uint64_t endIdx = parIdx+1 == smallerSrcCnt.size() ? smallerSrcCnt.size() : smallerSrcCnt[parIdx+1];
        for (auto i = startIdx; i < endIdx; i++) {
//            std::cerr << "s" << i << " ";
            if (i >= smallerSrc.size()) {
                std::cerr << "3 happened\n";
                std::exit(3);
            }
            auto child = smallerSrc[i] >> weightBits;
//            std::cerr << child << "\n";
            if (child >= visited.size()) {
                std::cerr << "4 happened\n";
                std::exit(3);
            }
            if (child >= parentbv.size()) {
                std::cerr << "5 happened\n";
                std::exit(3);
            }
            if (!visited[child]) {
                parentbv[child] = parIdx;
                dfs(child, parentbv, visited, activeLevel, remaining, level);
            }
        }
        if (parIdx >= greaterSrcCnt.size()) {
            std::cerr << "22 happened\n";
            std::exit(3);
        }
        startIdx = greaterSrcCnt[parIdx];
        endIdx = parIdx+1 == greaterSrcCnt.size() ? greaterSrcCnt.size() : greaterSrcCnt[parIdx+1];
        for (auto i = startIdx; i < endIdx; i++) {
//            std::cerr << "g" << i << " ";
            if (i >= greaterSrc.size()) {
                std::cerr << "6 happened\n";
                std::cerr << parIdx << " " << i << " " << endIdx << "\n";
                std::exit(3);
            }
            auto child = greaterSrc[i];
//            std::cerr << child << "\n";
            if (child >= visited.size()) {
                std::cerr << "7 happened\n";
                std::exit(3);
            }
            if (child >= parentbv.size()) {
                std::cerr << "8 happened\n";
                std::exit(3);
            }
            if (!visited[child]) {
                parentbv[child] = parIdx;
                dfs(child, parentbv, visited, activeLevel, remaining, level);
            }
        }
//        std::cerr << "\n";
    }

/*
    std::vector<colorIdType> dfsWalk(uint64_t idx, sdsl::bit_vector &visited) {
        visited[idx] = 1;
        uint64_t startIdx = smallerSrcCnt[idx];
        uint64_t endIdx = idx+1 == smallerSrcCnt.size() ? smallerSrcCnt.size() : smallerSrcCnt[idx+1];
        for (auto i = startIdx; i < endIdx; i++) {
            auto val = smallerSrc[i] >> weightBits;
            parentbv[i] = idx;
            dfsWalk()
            if (!visited[val]) {
                res.push_back(val);
            }
        }
        startIdx = greaterSrcCnt[idx];
        endIdx = idx+1 == greaterSrcCnt.size() ? greaterSrcCnt.size() : greaterSrcCnt[idx+1];
        for (auto i = startIdx; i < endIdx; i++) {
            if (!visited[greaterSrc[i]]) {
                res.push_back(greaterSrc[i]);
            }
        }
        return res;
    }
*/

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

    void kruskalMSF(AdjList * adjListPtr);

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

    void calcDeltasInParallel(uint32_t threadID, uint64_t deltaOffset,
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
