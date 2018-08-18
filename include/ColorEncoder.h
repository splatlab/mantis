//
// Created by Fatemeh Almodaresi on 8/17/18.
//

#ifndef MANTIS_COLORENCODER_H
#define MANTIS_COLORENCODER_H

#include <unordered_map>
#include "lru/lru.hpp"
#include "DeltaManager.h"

using LRUCacheMap =  LRU::Cache<uint64_t, std::vector<uint64_t>>;

class ColorEncoder {
public:
    ColorEncoder(uint64_t numSamplesIn) : numSamples(numSamplesIn),
                                          deltaM(numSamplesIn, 222500000, 8) {}
    void addEdge(uint64_t i, uint64_t j, uint32_t w);
    bool hasEdge(uint64_t i, uint64_t j);
    std::vector<uint64_t> buildColor(uint64_t eqid);
    DeltaManager deltaM;

private:
    uint64_t numSamples;

    struct pair_hash {
        template <class T1, class T2>
        std::size_t operator () (const std::pair<T1,T2> &p) const {
            auto h1 = std::hash<T1>{}(p.first);
            auto h2 = std::hash<T2>{}(p.second);

            // Mainly for demonstration purposes, i.e. works but is overly simple
            // In the real world, use sth. like boost.hash_combine
            return h1 ^ h2;
        }
    };
    std::unordered_map<std::pair<uint64_t, uint64_t>, uint32_t, pair_hash> edges;
    LRU::Cache <uint64_t, std::vector<uint64_t>> lru_cache;
    //DeltaManager deltas;
};
#endif //MANTIS_COLORENCODER_H
