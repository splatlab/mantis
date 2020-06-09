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
#include <parallel/algorithm>

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
#include "adjList.h"
#include "cqfMerger.h"

using LRUCacheMap =  LRU::Cache<uint64_t, std::vector<uint64_t>>;

//using LRUCacheMap = HPHP::ConcurrentScalableCache<uint64_t , std::vector<uint64_t >>;

using SpinLockT = std::mutex;

using FilterType = CQF<KeyObject>;

typedef uint64_t colorIdType;
typedef uint32_t weightType;
typedef unsigned __int128 uint128_t;


class TmpFileIterator {
    typedef std::pair<uint64_t , uint64_t > valType;
public:
    TmpFileIterator(std::string fileNameIn, uint64_t maxBufferSizeIn) :
    filename(fileNameIn), maxBufferSize(maxBufferSizeIn) {
        file.open(filename, std::ios::in | std::ios::binary);
        file.read(reinterpret_cast<char* >(&countOfItemsInFile), sizeof(countOfItemsInFile));
        next();
    }
    ~TmpFileIterator() {file.close();}
    TmpFileIterator(TmpFileIterator&& o) noexcept {
        filename = o.filename;
        fileIdx = o.fileIdx;
        vecIdx = o.vecIdx;
        countOfItemsInFile = o.countOfItemsInFile;
        maxBufferSize = o.maxBufferSize;
        buffer = o.buffer;
        o.buffer.clear();
        file.open(filename, std::ios::in | std::ios::binary);
        file.seekg(o.file.tellg());
        o.file.close();
    }

    bool next() {
        vecIdx++;
        if (end()) return false;
        if (vecIdx >= buffer.size()) {
            auto toFetch = std::min(maxBufferSize, countOfItemsInFile - fileIdx);
            fileIdx += toFetch;
            buffer.resize(toFetch);
            file.read(reinterpret_cast<char*>(buffer.data()), sizeof(valType)*toFetch);
            vecIdx = 0;
        }
        return true;
    }

    bool end() const {
        return fileIdx >= countOfItemsInFile and vecIdx >= buffer.size();
    }

    bool operator>(const TmpFileIterator &rhs) const {
        return get_val() > rhs.get_val();
    }

    const valType &get_val() const { return buffer[vecIdx]; }

    std::string get_filename() {return filename;}
private:
    std::string filename;
    std::ifstream file;
    uint64_t fileIdx{0};
    uint64_t vecIdx{0};
    uint64_t countOfItemsInFile{0};
    uint64_t maxBufferSize{0};
    std::vector<valType> buffer;
};

struct Minheap_edge {
    void push(TmpFileIterator *obj) {
        c.emplace_back(obj);
        std::push_heap(c.begin(), c.end(), [](auto &f, auto &s) {return (*f) > (*s);});
    }

    void pop() {
        std::pop_heap(c.begin(), c.end(), [](auto &f, auto &s) {return (*f) > (*s);});
        c.pop_back();
    }

    void replace_top(TmpFileIterator *obj) {
        c.emplace_back(obj);
        pop();
    }

    TmpFileIterator* top() {return c.front(); }

    bool empty() const { return c.empty(); }

private:
    std::vector<TmpFileIterator*> c;
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

// To represent Disjoint Sets
struct DisjointTrees {
    sdsl::int_vector<> els{};

    // Constructor.
    explicit DisjointTrees(uint64_t n) {
        // Allocate memory
//        this->n = n;
        els = sdsl::int_vector<>(n, 1, ceil(log2(n))+1); // 1 bit which is set if the node IS its own parent
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
    void merge(uint64_t x, uint64_t y) {
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

class MSTMerger {
public:
    MSTMerger(std::string prefix,
            spdlog::logger *logger,
            uint32_t numThreads,
            std::string prefix1,
            std::string prefix2);

    void mergeMSTs();

private:

    bool buildEdgeSets();

    bool calculateMSTBasedWeights();

    bool encodeColorClassUsingMST();

    void kruskalMSF(AdjList * adjListPtr);

    std::set<workItem> neighbors(CQF<KeyObject> &cqf, workItem n);

    bool exists(CQF<KeyObject> &cqf, dna::canonical_kmer e, uint64_t &eqid);

    uint32_t mstBasedHammingDist(uint64_t eqid1,
                                 uint64_t eqid2,
                                 MSTQuery *mst,
                                 LRUCacheMap &lru_cache,
                                 QueryStats &queryStats,
                                 std::unordered_map<uint64_t, std::vector<uint64_t>> &fixed_cache);

    void buildPairedColorIdEdgesInParallel(uint32_t threadId, CQF<KeyObject> &cqf,
                                           std::vector<std::pair<uint64_t, uint64_t>> &tmpEdges,
                                           uint64_t &curFileIdx,
                                           uint64_t &cnt, uint64_t &maxId, uint64_t &numOfKmers,
                                           spp::sparse_hash_map<std::pair<uint64_t , uint64_t >, uint64_t, Custom_Pair_Hasher > &popularEdges);


    void findNeighborEdges(CQF<KeyObject> &cqf, KeyObject &keyobj, std::vector<std::pair<uint64_t, uint64_t>> &edgeList,
                           spp::sparse_hash_map<std::pair<uint64_t , uint64_t >, uint64_t, Custom_Pair_Hasher > & popularEdges);

    void calcDeltasInParallel(uint32_t threadID, uint64_t &deltaOffset,
                              sdsl::int_vector<> &parentbv, sdsl::int_vector<> &deltabv,
                              sdsl::bit_vector &bbv,
                              bool isMSTBased,
                              AdjList * adjListPtr
                              );

    void buildMSTBasedColor(uint64_t eqid1, MSTQuery *mst1,
                            LRUCacheMap &lru_cache1, std::vector<uint64_t> &eq1,
                            QueryStats &queryStats,
                            std::unordered_map<uint64_t, std::vector<uint64_t>> &fixed_cache);

    std::vector<uint32_t> getMSTBasedDeltaList(uint64_t eqid1, uint64_t eqid2,
                                               MSTQuery * mstPtr,
                                               std::unordered_map<uint64_t, std::vector<uint64_t>>& fixed_cache,
                                               LRUCacheMap &lru_cache,
                                               QueryStats &queryStats,
                                               bool verbose=false);

    void edgePairSortUniq(std::vector<std::pair<uint64_t, uint64_t>> &edgeList) {
//        omp_set_dynamic(false);
//        omp_set_num_threads(nThreads);
        std::sort(std::execution::par_unseq,edgeList.begin(), edgeList.end(),
                             [](auto &e1, auto &e2) {
                                 return e1.first == e2.first ? e1.second < e2.second : e1.first < e2.first;
                             });
        edgeList.erase(std::unique(edgeList.begin(), edgeList.end(),
                                   [](auto &e1, auto &e2) {
                                       return e1.first == e2.first and e1.second == e2.second;
                                   }), edgeList.end());
    }
    inline uint64_t splitHarmonically(uint64_t size, uint64_t cnt, std::vector<uint64_t> &endIndices) {
        endIndices.resize(cnt);
        // alpha * (1 + 1/2 + ... + 1/cnt) <= alpha * (log(cnt)+1) = size
        // alpha = size/(log(cnt)+1)
        auto alpha = static_cast<uint64_t >(ceil(size / (log(cnt)+1)));
//        std::cerr << "alpha is " << alpha << " for size: " << size << " and count: " << cnt << "\n";
        size = 0;
        for (uint64_t i = 0; i < endIndices.size(); i++) {
            double ratio = 1.0/(i+1);
            auto val = static_cast<uint64_t >(ceil(alpha * ratio));
            size += val;
            endIndices[i] = size;
//            std::cerr  << startIndices[i] << " ";
        }
//        std::cerr << "\n new size: " << size << "\n";
        return size;
    }

    inline void removeIntermediateDiskFiles() {
        std::string sysCommand = "rm -r " + prefix + "tmp*";
        system(sysCommand.c_str());
        sysCommand = "rm -r " + prefix + "w*";
        system(sysCommand.c_str());
    }

    inline void storeMST(sdsl::int_vector<> &parentbv,
            sdsl::int_vector<> &deltabv,
            sdsl::bit_vector &bbv) {
        logger->info("Serializing data structures parentbv, deltabv, & bbv...");
        sdsl::store_to_file(parentbv, std::string(prefix + mantis::PARENTBV_FILE));
        sdsl::store_to_file(deltabv, std::string(prefix + mantis::DELTABV_FILE));
        sdsl::store_to_file(bbv, std::string(prefix + mantis::BOUNDARYBV_FILE));
        logger->info("Done Serializing.");
    }

    inline std::vector<uint64_t> findThreadWeightBoundaries(sdsl::int_vector<> &parentbv,
                                                            AdjList *adjListPtr) {
        std::vector<uint64_t> thread_deltaOffset_and_parentEnd(nThreads, 0);
        uint64_t idx{0};
        uint64_t bucketSize = std::ceil(parentbv.size() / (double)nThreads);
        uint64_t doubleCheckTotWeight{0};
        for (auto i = 1; i <= adjListPtr->smallerSrcStartIdx.size(); i++) {
            // Limit for adj edges of parent, are the start of the adj edges for next node.
            // So if parent = i-1, limit for adj edges of parent is startIdx[i]
            auto limit = i == adjListPtr->smallerSrcStartIdx.size() ?
                         adjListPtr->smallerSrcStartIdx.size() : adjListPtr->smallerSrcStartIdx[i];
            while (idx < limit) {
                auto par = static_cast<uint64_t>(i-1);
                auto val = adjListPtr->smallerSrc[idx];
                uint64_t weight = val & adjListPtr->weightMask;
                uint64_t child = val >> adjListPtr->weightBits;
                if (parentbv[par] == child) {
                    std::swap(par, child);
                } else if (parentbv[child] != par) {
                    std::cerr << "ERROR! Neither of the two nodes are the parent at index: " << idx << "\n" <<
                              "Expected: " << child << " <-> " << par << "\n"
                              << "Got: " << parentbv[par] << " for " << par <<
                              " and " << parentbv[child] << " for " << child << "\n";
                    std::exit(3);
                }
                thread_deltaOffset_and_parentEnd[std::min(static_cast<uint64_t >(nThreads-1), child / bucketSize)] += weight;
                doubleCheckTotWeight += weight;
                idx++;
            }
        }
        if (doubleCheckTotWeight != mstTotalWeight - 1) {
            std::cerr << "ERROR! Weights are not stored properly:\n" <<
                      "Expected: " << mstTotalWeight << " Got: " << doubleCheckTotWeight << "\n";
            std::exit(3);
        }
        for (auto i = 1; i < thread_deltaOffset_and_parentEnd.size(); i++) {
            thread_deltaOffset_and_parentEnd[i] += thread_deltaOffset_and_parentEnd[i-1];
        }
        for (auto i = thread_deltaOffset_and_parentEnd.size()-1; i > 0 ; i--) {
            thread_deltaOffset_and_parentEnd[i] = thread_deltaOffset_and_parentEnd[i-1];
//        std::cerr << "thr" << i << " " << thread_deltaOffset_and_parentEnd[i] << "\n";
        }
        thread_deltaOffset_and_parentEnd[0] = 0;
        return thread_deltaOffset_and_parentEnd;
    }

    std::string prefix;
    uint32_t numSamples = 0;
    uint32_t toBeMergedNumOfSamples[2];
    uint64_t k;
    uint64_t num_colorClasses = 0;
    uint64_t mstTotalWeight = 0;
    uint64_t zero = static_cast<colorIdType>(UINT64_MAX);
    std::vector<LRUCacheMap> lru_cache[2];//10000);
    std::unordered_map<uint64_t, std::vector<uint64_t>> fixed_cache[2];
    std::vector<QueryStats> queryStats[2];
    sdsl::int_vector<> colorPairs[2];
    sdsl::int_vector<> ccBits;
    std::vector<uint64_t> ccBitsBucketCnt;
    std::string prefixes[2];
    std::unique_ptr<MSTQuery> mst[2];
    spdlog::logger *logger{nullptr};
    uint32_t nThreads = 1;
    SpinLockT colorMutex, writeMutex;

    uint64_t numBlocks;
    uint64_t curFileIdx = 0;
    uint32_t maxWeightInFile{1000};
    uint64_t ccCnt[2];

};

#endif //MANTIS_MST_H
