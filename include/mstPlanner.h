//
// Created by Fatemeh Almodaresi on 2020-06-04.
//

#ifndef MANTIS_MSTPLANNER_H
#define MANTIS_MSTPLANNER_H


struct Cost {
    uint64_t numSteps{1};
    uint64_t numQueries{1};
};

class MSTPlanner {
public:
    MSTPlanner(std::string &prefix,
               MSTQuery *mst0,
               MSTQuery *mst1,
               sdsl::int_vector<> &colorPairs0,
               sdsl::int_vector<> &colorPairs1,
               uint64_t edgeFileCnt,
               uint64_t maxAllowedEdgesInBuffer,
               uint64_t nThreads,
               spdlog::logger *loggerIn,
               bool includingQueries = true) {
        logger = loggerIn;
        MSTQuery* mst[2];
        mst[0] = mst0;
        mst[1] = mst1;

        for (uint64_t mstIdx = 0; mstIdx < 2; mstIdx++) {
            mstCost[mstIdx].resize(mst[mstIdx]->parentbv.size());
            mstZero[mstIdx] = mst[mstIdx]->parentbv.size()-1;
        }
        // local step cost for each node is its degree (in + out)
        // query cost for each node is number of times it is called
        logger->info("Loading {} temp edge files.", edgeFileCnt);
        if (includingQueries) {
            uint64_t tmpedgecnt = 0;
            std::vector<TmpFileIterator> tis;
            Minheap_edge minheap;
            for (uint64_t i = 0; i < edgeFileCnt; i++) {
                tis.emplace_back(prefix + "tmp" + std::to_string(i), maxAllowedEdgesInBuffer);
            }
            for (uint64_t i = 0; i < edgeFileCnt; i++) {
                if (tis[i].end()) continue;
                minheap.push(&tis[i]);
            }
            uint64_t maxAllowedCnt = nThreads * 1000000;
            std::vector<std::pair<colorIdType, colorIdType>> uniqueEdges[2];
            uniqueEdges[0].reserve(maxAllowedCnt);
            uniqueEdges[1].reserve(maxAllowedCnt);
            auto val = minheap.top()->get_val();
            while (!minheap.empty()) {
                do {
                    auto cur = minheap.top();
                    val = cur->get_val();
                    if (cur->next()) {
                        minheap.replace_top(cur);
                    } else
                        minheap.pop();
                } while (!minheap.empty() && val == minheap.top()->get_val());
                if (val.first == val.second) {
                    logger->error("edge.first = edge.second while no self-edges are allowed: {}->{}", val.first, val.second);
                    std::exit(3);
                }
                std::pair<uint64_t , uint64_t> edge = std::make_pair(colorPairs0[val.first], colorPairs0[val.second]);
                if (edge.first != edge.second)
                    uniqueEdges[0].push_back(edge);
                if (edge.first > mstZero[0] or edge.second > mstZero[0]) {
                    std::cerr << "FFFFUUUCK0 " <<
                    mstZero[0] << " " << edge.first << " " << edge.second
                            << " " << val.first << " " << val.second << "\n";
                    std::exit(3);
                }
                edge = std::make_pair(colorPairs1[val.first], colorPairs1[val.second]);
                if (edge.first > mstZero[1] or edge.second > mstZero[1]) {
                    std::cerr << "FFFFUUUCK1 " << mstZero[1] << " " << edge.first << " " << edge.second << "\n";
                    std::exit(3);
                }
                if (edge.first != edge.second)
                    uniqueEdges[1].push_back(edge);
                if (uniqueEdges[0].size() > maxAllowedCnt or minheap.empty()) {
                    for (uint64_t mstIdx = 0; mstIdx < 2; mstIdx++) {
//                        omp_set_dynamic(false);
//                        omp_set_num_threads(nThreads);
                        std::sort(
                                uniqueEdges[mstIdx].begin(), uniqueEdges[mstIdx].end(),
                                             [](auto &e1, auto &e2) {
                                                 return e1.first == e2.first ? e1.second < e2.second : e1.first <
                                                                                                       e2.first;
                                             });
                        uniqueEdges[mstIdx].erase(
                                std::unique(uniqueEdges[mstIdx].begin(), uniqueEdges[mstIdx].end(),
                                                              [](auto &e1, auto &e2) {
                                                                  return e1.first == e2.first and
                                                                         e1.second == e2.second;
                                                              }), uniqueEdges[mstIdx].end());
                        for (auto &e: uniqueEdges[mstIdx]) {
                            if (e.first != mstZero[mstIdx]) {
                                mstCost[mstIdx][e.first].numQueries++;
                            } else if (e.second != mstZero[mstIdx]) {
                                mstCost[mstIdx][e.second].numQueries++;
                            }
                        }
                        uniqueEdges[mstIdx].clear();
                    }

                }
            }
        }
        logger->info("Done setting the local costs");
    }

    bool plan(uint64_t mstIdx, MSTQuery *mst, std::unordered_map<uint64_t, std::vector<uint64_t>> &fixed_cache) {
        logger->info("Planner for mantis {}", mstIdx);
        LRUCacheMap lru_cache(1000);
        nonstd::optional<uint64_t> dummy{nonstd::nullopt};
        QueryStats dummyStats;
        std::vector<colorIdType> colorsInCache;
        std::vector<std::vector<colorIdType>> children(mst->parentbv.size());

        for (uint64_t i = 0; i < mstZero[mstIdx]; i++) {
            if (mst->parentbv[i] >= children.size()) {
                logger->error("A parent index should be in the range of valid indices. "
                              "parentIdx of {}:{}, maxIdx:{}",mstIdx,  mst->parentbv[i], children.size());
                std::exit(3);
            }
            children[mst->parentbv[i]].push_back(i);
        }
        // recursive planner
        uint64_t cntr{0};
        planRecursively(mstZero[mstIdx], mstCost[mstIdx], children, colorsInCache, cntr);
        logger->info("fixed cache size is : {}", colorsInCache.size());
        // fillout fixed_cache[0]
        // walking in reverse order on colorsInCache is intentional
        // because as we've filled out the colorsInCache vector in a post-order traversal of MST,
        // if we walk in from end to the beginning we can always guarantee that we've already put the ancestors
        // for the new colors we want to put in the fixed_cache and that will make the whole color bv construction
        // FASTER!!
        for (int64_t idx = colorsInCache.size() - 1; idx >= 0; idx--) {
            auto setbits = mst->buildColor(colorsInCache[idx], dummyStats, &lru_cache, nullptr, &fixed_cache, dummy);
            fixed_cache[colorsInCache[idx]] = setbits;
        }
    }
private:

    void  planRecursively(uint64_t nodeId,
                                    std::vector<Cost> &mstCostv,
                                    std::vector<std::vector<colorIdType>> &children,
                                    std::vector<colorIdType> &colorsInCache,
                                    uint64_t &cntr) {

        for (auto c : children[nodeId]) {
            planRecursively(c, mstCostv, children, colorsInCache, cntr);
        }

        std::unordered_set<colorIdType> localCache;
        uint64_t numQueries{0}, numSteps{0};
        colorIdType childWithMaxAvg{0};
        float maxChildAvg{0};
        do {
//        auto tmp = numQueries == 0 ? 0 : (float)numSteps/(float)numQueries;
//        std::cerr << "\rnodeId " << nodeId << " " << tmp << " " << numSteps << " " << numQueries << " " <<localCache.size();
            if (maxChildAvg != 0) {
                localCache.insert(childWithMaxAvg);
//                std::cerr << childWithMaxAvg << ":" << maxChildAvg << "\n";
                childWithMaxAvg = 0;
                maxChildAvg = 0;
            }

            numQueries = mstCostv[nodeId].numQueries;
            numSteps = mstCostv[nodeId].numSteps;
            for (auto c : children[nodeId]) {
                if (localCache.find(c) == localCache.end()) {
                    numQueries += mstCostv[c].numQueries;
                    // steps in parent += steps in children + 1 for each child
                    numSteps += mstCostv[c].numSteps + mstCostv[c].numQueries;
                    float currChildAvg =
                            mstCostv[c].numQueries == 0 ? 0 : (float) mstCostv[c].numSteps / (float) mstCostv[c].numQueries;
                    if (mstCostv[c].numSteps != 0 and
                        (currChildAvg > maxChildAvg or
                         (currChildAvg == maxChildAvg and mstCostv[childWithMaxAvg].numQueries < mstCostv[c].numQueries))) {
                        maxChildAvg = currChildAvg;
                        childWithMaxAvg = c;
                    }
                }
            }
        } while (numQueries != 0 and (float) numSteps / (float) numQueries > avgCostThreshold);
        mstCostv[nodeId].numQueries = numQueries;
        mstCostv[nodeId].numSteps = numSteps;
        for (auto c : localCache) {
            colorsInCache.push_back(c);
        }
        if (++cntr % 1000000 == 0)
            std::cerr << "\r" << cntr++;
    }



    std::vector<Cost> mstCost[2];
    uint64_t mstZero[2]{};
    spdlog::logger *logger;
    const uint64_t avgCostThreshold = 16;
};
#endif //MANTIS_MSTPLANNER_H
