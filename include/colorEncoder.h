//
// Created by Fatemeh Almodaresi on 8/17/18.
//

#ifndef MANTIS_COLORENCODER_H
#define MANTIS_COLORENCODER_H

#include <unordered_map>
#include <unordered_set>
#include "lru/lru.hpp"
#include "deltaManager.h"
#include "sdsl/bit_vectors.hpp"
#include "cqf.h"
#include "canonKmer.h"
#include "hashutil.h"

using LRUCacheMap =  LRU::Cache<uint64_t, std::vector<uint64_t>>;
typedef std::pair<duplicated_dna::canonical_kmer, uint64_t> node;
constexpr uint64_t zero = 0;

struct Edge {
    uint64_t parent;
    uint64_t child;
    uint64_t weight;

    Edge() {
        parent = 0; child = 0; weight = 0;
    }
    Edge(uint64_t inParent, uint64_t inChild, uint64_t inWeight) :
            parent(inParent), child(inChild), weight(inWeight) {}
};

struct pair_hash {
    template<class T1, class T2>
    std::size_t operator()(const std::pair<T1, T2> &p) const {
        auto h1 = std::hash<T1>{}(p.first);
        auto h2 = std::hash<T2>{}(p.second);

        // Mainly for demonstration purposes, i.e. works but is overly simple
        // In the real world, use sth. like boost.hash_combine
        return h1 ^ h2;
    }
};

class ColorEncoder {
public:

    ColorEncoder(uint64_t numSamplesIn,
                 CQF<KeyObject> &cqfIn,
                 uint64_t approximateClrClsesIn,
                 uint64_t approximateDeltaCntPerClrCls = 8) :
            numSamples(numSamplesIn),
            cqf(cqfIn),
            bvSize(approximateClrClsesIn),
            parentbv(bvSize, 0, log2((double)bvSize)+5),//TODO take care of this constant!!
            deltaM(numSamplesIn, bvSize, approximateDeltaCntPerClrCls),
            colorClsCnt(1), // start with the dummy node
            lru_cache(100000)
            {
                std::cerr << "\nColorEncoder Constructor:  bvSize: "
                             << bvSize << " parent size: " << parentbv.size()
                          << " colorClsCnt: " << colorClsCnt << "\n";
            }


    bool addColorClass(uint64_t kmer, uint64_t eqId, const sdsl::bit_vector &bv);

    bool serialize(std::string prefix);

private:
    uint64_t numSamples;
    uint64_t bvSize;
    uint64_t colorClsCnt;
    sdsl::int_vector<> parentbv;
    DeltaManager deltaM;
    CQF<KeyObject> &cqf;
    int k = 20;
    std::unordered_map<std::pair<uint64_t, uint64_t>, uint32_t, pair_hash> edges;
    LRUCacheMap lru_cache;

    std::vector<uint64_t> buildColor(uint64_t eqid);

    std::vector<uint64_t> buildColor(const sdsl::bit_vector &bv);

    bool updateMST(uint64_t n1, uint64_t n2, std::vector<uint64_t> deltas);

    std::vector<uint64_t> hammingDist(uint64_t i, uint64_t j);

    std::unordered_set<uint64_t> neighbors(duplicated_dna::canonical_kmer n);

    bool exists(duplicated_dna::canonical_kmer e, uint64_t &eqid);

    std::pair<Edge, Edge> maxWeightsTillLCA(uint64_t n1, uint64_t n2);

    void addEdge(uint64_t i, uint64_t j, uint32_t w);

    bool hasEdge(uint64_t i, uint64_t j);

    uint32_t getEdge(uint64_t i, uint64_t j);
};

#endif //MANTIS_COLORENCODER_H
