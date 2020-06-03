//
// Created by Fatemeh Almodaresi on 2020-05-13.
//

#include <cmath>
#include <mantisconfig.hpp>
#include "adjList.h"

AdjList::AdjList(std::string prefixIn, uint64_t numColorClasses, uint64_t numSamples) : prefix(prefixIn), adjListFile(prefix+mantis::TEMP_MST_ADJ_FILE) {
    weightBits = static_cast<uint64_t>(ceil(log2(numSamples)));
    weightMask = (1ULL << weightBits) - 1;
    smallerSrc = sdsl::int_vector<>(numColorClasses, 0, ceil(log2(numColorClasses))+weightBits);
    smallerSrcStartIdx = sdsl::int_vector<>(numColorClasses, 0, ceil(log2(numColorClasses)));
    greaterSrc = sdsl::int_vector<>(numColorClasses, 0, ceil(log2(numColorClasses)));
    greaterSrcStartIdx = sdsl::int_vector<>(numColorClasses, 0, ceil(log2(numColorClasses)));
}

/**
     * IMportant NOTE:: Assumes src < dest!!
     * @param src smaller vertex of the edge
     * @param dest greater vertex of the edge
     * @param weight of the edge
     */
void AdjList::storeEdge(uint64_t src, uint64_t dest, uint64_t weight) {
    if (src > dest) {
        std::cerr << "WARNING!! Expect src < dest but " << src << " > " << dest << "\n";
        std::swap(src, dest);
    }
    adjListFile << src << " " << dest << " " << weight << "\n";
    smallerSrcStartIdx[src] = smallerSrcStartIdx[src]+1;
    greaterSrcStartIdx[dest] = greaterSrcStartIdx[dest]+1;
//    std::cerr << "edge " << idx << " " << src << " " << dest << " " << smallerSrcStartIdx[src] << " " << greaterSrcStartIdx[dest] << "\n";
}

void AdjList::loadCompactedAdjList() {
    for (auto i = 1; i < smallerSrcStartIdx.size(); i++) {
        smallerSrcStartIdx[i] = smallerSrcStartIdx[i-1] + smallerSrcStartIdx[i];
        greaterSrcStartIdx[i] = greaterSrcStartIdx[i-1] + greaterSrcStartIdx[i];
    }

    adjListFile.close();
    std::ifstream adjListIn(prefix+mantis::TEMP_MST_ADJ_FILE);
    uint64_t src, dest, weight;
    adjListIn >> src >> dest >> weight;
    while (adjListIn.good()) {
        smallerSrcStartIdx[src] = smallerSrcStartIdx[src] - 1;
        smallerSrc[smallerSrcStartIdx[src]] = ((uint128_t)dest << weightBits) | weight;

        greaterSrcStartIdx[dest] = greaterSrcStartIdx[dest] - 1;
        greaterSrc[greaterSrcStartIdx[dest]] = src;
        adjListIn >> src >> dest >> weight;
    }
    adjListIn.close();
    // system("rm " + prefix + mantis::TEMP_MST_ADJ_FILE);
}

/**
 * dfs + bfs
 */
void AdjList::hybridTreeWalk(uint64_t root, sdsl::int_vector<> &parentbv, sdsl::bit_vector &visited, uint64_t level) {
    sdsl::bit_vector currentLevel(parentbv.size(), 0);
    uint64_t remaining = 0;
    parentbv[root] = root; // and it's its own parent (has no parent)
    uint64_t parIdx = root;
    boundedDfs(parIdx, parentbv, visited, currentLevel, remaining);
//    std::cerr << "remaining " << remaining << "\n";
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
//                        std::cerr << "in if remaining " << remaining << " ";
                        boundedDfs(parIdx, parentbv, visited, currentLevel, remaining);
//                        std::cerr << " after remaining " << remaining << "\n";
                    }
                }
            }
            auto newWrd = currentLevel.get_int(idx, wrdLen);
            currentLevel.set_int(idx, newWrd ^ wrd, wrdLen); // this is so cool and smart :D
            idx += wrdLen;
        }
//        std::cerr << "remaining " << remaining << "\n";
    }
}
/**
 * DFS walk instead of BFS
 */
void AdjList::boundedDfs(uint64_t parIdx,
                sdsl::int_vector<> &parentbv,
                sdsl::bit_vector &visited,
                sdsl::bit_vector &activeLevel,
                uint64_t &remaining,
                uint64_t level) {
    if (level == 200) {
        activeLevel[parIdx] = 1;
        remaining++;
//        std::cerr << "200 remaining " << remaining << "\n";
        return;
    }
    if (parIdx >= visited.size()) {
        std::cerr << "1 happened\n";
        std::exit(3);
    }
    if (visited[parIdx] == 0) {
        visitedCnt++;
    }
//    std::cerr << "l" << level << " => " << parIdx << " ";
    visited[parIdx] = 1;
    level++;

    if (parIdx >= smallerSrcStartIdx.size()) {
        std::cerr << "2 happened\n";
        std::exit(3);
    }
    uint64_t startIdx = smallerSrcStartIdx[parIdx];
    uint64_t endIdx = parIdx+1 == smallerSrcStartIdx.size() ? smallerSrc.size() : smallerSrcStartIdx[parIdx+1];
    auto cnt = endIdx - startIdx;
    auto tmp1 = greaterSrcStartIdx[parIdx];
    auto tmp2 = parIdx+1 == greaterSrcStartIdx.size() ? greaterSrcStartIdx.size() : greaterSrcStartIdx[parIdx+1];
//    std::cerr << " child cnt: " << cnt + (tmp2 - tmp1) - 1 << "\n";
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
            boundedDfs(child, parentbv, visited, activeLevel, remaining, level);
        }
    }
    if (parIdx >= greaterSrcStartIdx.size()) {
        std::cerr << "22 happened\n";
        std::exit(3);
    }
    startIdx = greaterSrcStartIdx[parIdx];
    endIdx = parIdx+1 == greaterSrcStartIdx.size() ? greaterSrc.size() : greaterSrcStartIdx[parIdx+1];
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
            boundedDfs(child, parentbv, visited, activeLevel, remaining, level);
        }
    }
//        std::cerr << "\n";
}
