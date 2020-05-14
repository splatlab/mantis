//
// Created by Fatemeh Almodaresi on 2020-05-13.
//

#ifndef MANTIS_ADJLIST_H
#define MANTIS_ADJLIST_H


#include <cstdint>
#include <iostream>
#include <memory>
#include <fstream>
#include "sdsl/bit_vectors.hpp"

typedef unsigned __int128 uint128_t;

/**
 * Adjacency List
 * Weight is stored within the smallerSrc vector
 */
class AdjList {
private:
    std::string prefix;
    std::ofstream adjListFile;
    uint64_t numCCs{0}, numSamples{0};
public:
    // weight is in the smallerSrc
    sdsl::int_vector<> smallerSrc;
    sdsl::int_vector<> smallerSrcCnt;
    sdsl::int_vector<> greaterSrc;
    sdsl::int_vector<> greaterSrcCnt;

    uint64_t idx{0}, visitedCnt{0};
    uint64_t weightBits;
    uint64_t weightMask;

    AdjList(std::string prefix, uint64_t numColorClasses, uint64_t numSamples);

    void storeEdge(uint64_t src, uint64_t dest, uint64_t weight);

    void loadCompactedAdjList();

    /**
     * dfs + bfs
     */
    void hybridTreeWalk(uint64_t root, sdsl::int_vector<> &parentbv, sdsl::bit_vector &visited, uint64_t level=0);

    void boundedDfs(uint64_t parIdx,
                    sdsl::int_vector<> &parentbv,
                    sdsl::bit_vector &visited,
                    sdsl::bit_vector &activeLevel,
                    uint64_t &remaining,
                    uint64_t level = 0);
};


#endif //MANTIS_ADJLIST_H
