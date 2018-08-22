//
// Created by Fatemeh Almodaresi on 8/17/18.
//

#ifndef MANTIS_COLORENCODER_H
#define MANTIS_COLORENCODER_H

#include <unordered_map>
#include <unordered_set>
#include "lru/lru.hpp"
#include "deltaManager.h"
#include "bitvector.h"
#include "cqf.h"
#include "canonKmer.h"

using LRUCacheMap =  LRU::Cache<uint64_t, std::vector<uint64_t>>;
typedef std::pair<dna::canonical_kmer, uint64_t> node;
constexpr uint64_t zero = 0;

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
            colorClsCnt(1) // start with the dummy node
            {}


    bool addColorClass(uint64_t kmer, uint64_t eqId, const std::vector<uint32_t> bv);

private:
    uint64_t numSamples;
    uint64_t bvSize;
    uint64_t colorClsCnt;
    sdsl::int_vector<> parentbv;
    DeltaManager deltaM;
    CQF<KeyObject> &cqf;
    int k = 20;
    std::unordered_map<std::pair<uint64_t, uint64_t>, uint32_t, pair_hash> edges;
    LRU::Cache<uint64_t, std::vector<uint64_t>> lru_cache;

    void addEdge(uint64_t i, uint64_t j, uint32_t w);

    bool hasEdge(uint64_t i, uint64_t j);

    std::vector<uint64_t> buildColor(uint64_t eqid);

    bool updateMST(uint64_t n1, uint64_t n2, std::vector<uint64_t> deltas);

    std::vector<uint64_t> hammingDist(uint64_t i, uint64_t j);

    std::unordered_set<uint64_t> neighbors(dna::canonical_kmer n);

    bool exists(dna::canonical_kmer e, uint64_t &eqid);

    void maxWeightsTillLCA(uint64_t n1, uint64_t n2,
                           uint64_t &w1, uint64_t &w2,
                           uint64_t &p1, uint64_t &p2);
};

#endif //MANTIS_COLORENCODER_H
