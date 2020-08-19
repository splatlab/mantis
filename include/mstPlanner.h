//
// Created by Fatemeh Almodaresi on 2020-06-04.
//

#ifndef MANTIS_MSTPLANNER_H
#define MANTIS_MSTPLANNER_H

#include "tbb/parallel_sort.h"
static constexpr uint64_t MAX_ALLOWED_TMP_EDGES_IN_FILE{64000000};//{130000000}; // Up to 2G


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
        if (includingQueries) {
            logger->info("Loading {} temp edge files.", edgeFileCnt);
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
            uint64_t maxAllowedCnt = MAX_ALLOWED_TMP_EDGES_IN_FILE;
            std::vector<std::pair<colorIdType, colorIdType>> uniqueEdges[2];
            uniqueEdges[0].reserve(maxAllowedCnt);
            uniqueEdges[1].reserve(maxAllowedCnt);
            logger->info("Walking all the edges to find popular nodes.");
            auto val = minheap.top()->get_val();
            if (val.first >= colorPairs0.size() or val.second >= colorPairs0.size()) {
                logger->error("Reading a wrong edge from the tmp file: {}->{} while max nodeId:{}", val.first, val.second, colorPairs0.size());
//                std::cerr << minheap.top()->get_filename() << " " << minheap.top()->countOfItemsInFile << " " << minheap.top()->fileIdx << " " << minheap.top()->vecIdx << " " << minheap.top()->buffer.size() << "\n";
                std::exit(3);
            }
            while (!minheap.empty()) {
                do {
                    auto cur = minheap.top();
                    val = cur->get_val();
                    if (val.first >= colorPairs0.size() or val.second >= colorPairs0.size()) {
                        logger->error("Reading a wrong edge from the tmp file: {}->{} while max nodeId:{}", val.first, val.second, colorPairs0.size());
//                        std::cerr << minheap.top()->get_filename() << " " << minheap.top()->countOfItemsInFile << " " << minheap.top()->fileIdx << " " << minheap.top()->vecIdx << " " << minheap.top()->buffer.size() << "\n";
                        std::exit(3);
                    }
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
                            << " colorPair0.size():" << colorPairs0.size() << " " << val.first << " " << val.second << "\n";
                    std::exit(3);
                }
                edge = std::make_pair(colorPairs1[val.first], colorPairs1[val.second]);
                if (edge.first > mstZero[1] or edge.second > mstZero[1]) {
                    std::cerr << "FFFFUUUCK1 "
                              << mstZero[1] << " " << edge.first << " " << edge.second
                              << " colorPair1.size():" << colorPairs1.size() << " " << val.first << " " << val.second << "\n";
                    std::exit(3);
                }
                if (edge.first != edge.second)
                    uniqueEdges[1].push_back(edge);
                if (uniqueEdges[0].size() > maxAllowedCnt or minheap.empty()) {
                    for (uint64_t mstIdx = 0; mstIdx < 2; mstIdx++) {
//                        omp_set_dynamic(false);
//                        omp_set_num_threads(nThreads);
                        tbb::parallel_sort( uniqueEdges[mstIdx].begin(), uniqueEdges[mstIdx].end(),
                                            [](auto &e1, auto &e2) {
                                                return e1.first == e2.first ? e1.second < e2.second : e1.first <
                                                                                                      e2.first;
                                            });
//                        std::sort(std::execution::par_unseq,
//                                uniqueEdges[mstIdx].begin(), uniqueEdges[mstIdx].end(),
//                                             [](auto &e1, auto &e2) {
//                                                 return e1.first == e2.first ? e1.second < e2.second : e1.first <
//                                                                                                       e2.first;
//                                             });
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

    void constructStaticCache(uint64_t mstIdx, MSTQuery *mst, std::unordered_map<uint64_t, std::vector<uint64_t>> &fixed_cache) {
        std::vector<colorIdType> colorsInCache;
        plan(mstIdx, mst, colorsInCache);
        fillCache(mst, colorsInCache, fixed_cache);
    }

private:

    void  planRecursively(uint64_t nodeId,
                                    std::vector<Cost> &mstCostv,
                                    std::vector<uint64_t> &childStartIdx,
                                    std::vector<uint64_t> &children,
                                    std::vector<colorIdType> &colorsInCache,
                                    uint64_t &cntr) {

        for (auto startIdx = childStartIdx[nodeId]; startIdx < childStartIdx[nodeId+1]; startIdx++) {
            auto c = children[startIdx];
            planRecursively(c, mstCostv, childStartIdx, children, colorsInCache, cntr);
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
            for (auto startIdx = childStartIdx[nodeId]; startIdx < childStartIdx[nodeId+1]; startIdx++) {
                auto c = children[startIdx];
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
            std::cerr << "\r" << cntr++ << " nodes processed.";
    }

    void plan(uint64_t mstIdx, MSTQuery *mst, std::vector<colorIdType> &colorsInCache) {
        logger->info("Planner for mantis {}", mstIdx);
        std::vector<uint64_t> childStartIdx(mst->parentbv.size()+1, 0);
        std::vector<uint64_t> children(mst->parentbv.size(), 0);
//        std::vector<std::vector<colorIdType>> children(mst->parentbv.size());

        for (uint64_t i = 0; i < mstZero[mstIdx]; i++) {
            if (mst->parentbv[i] >= childStartIdx.size()) {
                logger->error("A parent index should be in the range of valid indices. "
                              "parentIdx of {}:{}, maxIdx:{}",mstIdx,  mst->parentbv[i], childStartIdx.size());
                std::exit(3);
            }
            childStartIdx[mst->parentbv[i]]++;
        }
        for (uint64_t i = 1; i < childStartIdx.size(); i++) {
            childStartIdx[i]+=childStartIdx[i-1];
        }
        if (childStartIdx.back() != mst->parentbv.size()-1) {
            logger->error("Total#ofChildren != Total#ofNodes-1 in MST {} ( {} != {} ). Is the MST fully connected?", mstIdx, childStartIdx.back(), mst->parentbv.size()-1);
            std::exit(3);
        }
        logger->info("Done setting startIndices for parent->children vector.");
        for (uint64_t i = 0; i < mstZero[mstIdx]; i++) {
            if (mst->parentbv[i] >= children.size()) {
                logger->error("A parent index should be in the range of valid indices. "
                              "parentIdx of {}:{}, maxIdx:{}",mstIdx,  mst->parentbv[i], children.size());
                std::exit(3);
            }
            uint64_t parentIdx = mst->parentbv[i];
            childStartIdx[parentIdx]--;
            if (childStartIdx[parentIdx] >= children.size()) {
                logger->error("Out of range index for the parent->children vector "
                              "for MST {} => startIndex:{}, vectorSize:{}",mstIdx,  childStartIdx[parentIdx], children.size());
                std::exit(3);
            }
            uint64_t startIdx = childStartIdx[parentIdx];
            children[startIdx] = i;
        }
        logger->info("Done constructing the parent->children vector.");
        // recursive planner
        uint64_t cntr{0};
        planRecursively(mstZero[mstIdx], mstCost[mstIdx], childStartIdx, children, colorsInCache, cntr);
        childStartIdx.clear();
        childStartIdx.shrink_to_fit();
        children.clear();
        children.shrink_to_fit();
        mstCost[mstIdx].clear();
        mstCost[mstIdx].shrink_to_fit();
        std::cerr << "\r";
        logger->info("Fixed cache size is : {}", colorsInCache.size());
    }

    void fillCache(MSTQuery *mst, std::vector<colorIdType> &colorsInCache, std::unordered_map<uint64_t, std::vector<uint64_t>> &fixed_cache) {
        LRUCacheMap lru_cache(1000);
        nonstd::optional<uint64_t> dummy{nonstd::nullopt};
        QueryStats dummyStats;
        // fillout fixed_cache[0]
        // walking in reverse order on colorsInCache is intentional
        // because as we've filled out the colorsInCache vector in a post-order traversal of MST,
        // if we walk in from end to the beginning we can always guarantee that we've already put the ancestors
        // for the new colors we want to put in the fixed_cache and that will make the whole color bv construction
        // FASTER!!
        logger->info("Filling the cache for planned color class IDs..");
        for (int64_t idx = colorsInCache.size() - 1; idx >= 0; idx--) {
            auto setbits = mst->buildColor(colorsInCache[idx], dummyStats, &lru_cache, nullptr, &fixed_cache, dummy);
            fixed_cache[colorsInCache[idx]] = setbits;
        }
        logger->info("Done filling the planned color class cache.");
    }

    std::vector<Cost> mstCost[2];
    uint64_t mstZero[2]{};
    spdlog::logger *logger;
    const uint64_t avgCostThreshold = 16;
};
#endif //MANTIS_MSTPLANNER_H
