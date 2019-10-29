#include <utility>

//
// Created by Fatemeh Almodaresi.
//
#include <string>
#include <sstream>
#include <cstdio>
#include <stdio.h>
#include <unordered_map>
#include <iterator>
#include <algorithm>

#include "MantisFS.h"
#include "mstMerger.h"
#include "ProgOpts.h"

#define MAX_ALLOWED_TMP_EDGES 31250000
#define LOWBIT_MASK 0xFFFFFFFF00000000
#define HIGHBIT_MASK 0x00000000FFFFFFFF
#define MAX_ALLOWED_BLOCK_SIZE 100
#define SHIFTBITS 32

static constexpr uint64_t NUMBlocks = 1 << 16;

MSTMerger::MSTMerger(CQF<KeyObject> *cqfIn, std::string prefixIn, spdlog::logger *loggerIn, uint32_t numThreads,
                     std::string prefixIn1, std::string prefixIn2, uint64_t numColorBuffersIn) :
        /*cqf(cqfIn),*/ prefix(std::move(prefixIn)),
        prefix1(std::move(prefixIn1)), prefix2(std::move(prefixIn2)),
        nThreads(numThreads),
        num_of_ccBuffers(numColorBuffersIn) {
    eqclass_files.resize(num_of_ccBuffers);
    logger = loggerIn;//.get();

    std::string cqf_file = std::string(prefix + mantis::CQF_FILE);
    cqf = new CQF<KeyObject>(cqf_file, CQF_FREAD);

    // Make sure the prefix is a full folder
    if (prefix.back() != '/') {
        prefix.push_back('/');
    }
    if (prefix1.back() != '/') {
        prefix1.push_back('/');
    }
    if (!mantis::fs::DirExists(prefix1.c_str())) {
        logger->error("Index parent directory for first mst, {}, does not exist", prefix1);
        std::exit(1);
    }
    if (prefix2.back() != '/') {
        prefix2.push_back('/');
    }
    if (!mantis::fs::DirExists(prefix2.c_str())) {
        logger->error("Index parent directory for second mst, {}, does not exist", prefix2);
        std::exit(1);
    }
    // make the output directory if it doesn't exist
    /*if (!mantis::fs::DirExists(prefix.c_str())) {
        logger->error("Index parent directory {} does not exist", prefix);
        std::exit(1);
    }
*/

    queryStats1.resize(nThreads);
    queryStats2.resize(nThreads);
    for (uint64_t t = 0; t < nThreads; t++) {
        lru_cache1.emplace_back(1000);
        lru_cache2.emplace_back(1000);
    }

    std::string sample_file = prefix1 + mantis::SAMPLEID_FILE;//(prefix.c_str() , mantis::SAMPLEID_FILE);
    std::ifstream sampleid(sample_file);
    std::string tmp;
    while (sampleid >> tmp >> tmp) {
        numSamples++;
    }
    sampleid.close();
    numOfFirstMantisSamples = numSamples; // this line is important!

    sample_file = prefix2 + mantis::SAMPLEID_FILE;//(prefix.c_str() , mantis::SAMPLEID_FILE);
    sampleid.open(sample_file);
    while (sampleid >> tmp >> tmp) {
        numSamples++;
    }
    secondMantisSamples = numSamples - numOfFirstMantisSamples;
    sampleid.close();

    logger->info("# of experiments: {}", numSamples);
    logger->info("{}, {}, {}, {}", numThreads, num_of_ccBuffers, num_colorClasses, num_edges);
}

/**
 * Builds an MST consists of 3 main steps:
 * 1. construct the color graph for all the colorIds derived from dbg
 *      This phase just requires loading the CQF
 * 2. calculate the weights of edges in the color graph
 *      This phase requires at most two buffers of color classes
 * 3. find MST of the weighted color graph
 */
void MSTMerger::mergeMSTs() {
    buildEdgeSets();
    calculateMSTBasedWeights();
    encodeColorClassUsingMST();
}

/**
 * iterates over all elements of CQF,
 * find all the existing neighbors, and build a color graph based on that
 * @return true if the color graph build was successful
 */
bool MSTMerger::buildEdgeSets() {
    std::vector<spp::sparse_hash_set<Edge, edge_hash>> edgesetList;
    edgeBucketList.resize(num_of_ccBuffers * num_of_ccBuffers);
    edgesetList.resize(num_of_ccBuffers * num_of_ccBuffers);

    logger->info("Reading colored dbg from disk.");
    std::string cqf_file(prefix + mantis::CQF_FILE);
    k = cqf->keybits() / 2;
    logger->info("Done loading cdbg. k is {}", k);
    logger->info("Iterating over cqf & building edgeSet ...");
    // max possible value and divisible by 64
    uint64_t maxId{0}, numOfKmers{0};

    // build color class edges in a multi-threaded manner
    std::vector<std::ofstream> blockFiles(NUMBlocks);
    std::vector<uint64_t> blockCnt(NUMBlocks, 0);
    blockMutex.resize(NUMBlocks);
    for (auto &b: blockMutex) { b = new std::mutex(); }
    std::queue<uint64_t> blockIds;
    for (uint64_t blockId = 0; blockId < NUMBlocks; blockId++) {
        blockFiles[blockId].open(prefix + "/tmpblock" + std::to_string(blockId),
                                 std::ios::out | std::ios::binary);
//        blockFiles[blockId].write(reinterpret_cast<char *>(&blockCnt[blockId]), sizeof(blockCnt[blockId]));
        blockIds.push(blockId);
    }
    std::vector<std::thread> threads;
    for (uint32_t i = 0; i < nThreads; ++i) {
        threads.emplace_back(std::thread(&MSTMerger::writePotentialColorIdEdgesInParallel, this, i,
                                         std::ref(*cqf), std::ref(blockIds), std::ref(blockFiles), std::ref(blockCnt)));
    }
    for (auto &t : threads) { t.join(); }
    logger->info("Done with writing down the block files.");
    threads.clear();
    uint64_t removed{0};
    std::vector<bool> blockValid(NUMBlocks, true);
    for (auto blockId = 0; blockId < NUMBlocks; blockId++) {
//        blockFiles[blockId].close();
        if (blockCnt[blockId] == 0) {
//            blockFiles[blockId].close();
            std::string file = prefix + "/tmpblock" + std::to_string(blockId);
            std::remove(file.c_str());
            blockValid[blockId] = false;
            removed++;
        }
        if (blockId % 1000 == 0) {
            std::cerr << "\r" << "out of " << blockId << " files " << removed << " were removed.";
        }
    }
//    blockFiles.clear();
    std::cerr << "\r";
    logger->info("{} block files were empty and removed.", removed);

    std::vector<std::ifstream> readBlockFiles(NUMBlocks);
    // blockIds should be empty at this point, so we refill it
    for (uint64_t blockId = 0; blockId < NUMBlocks; blockId++) {
        std::string blockFileName = prefix + "/tmpblock" + std::to_string(blockId);
//        if (mantis::fs::FileExists(blockFileName.c_str())) {
        if (blockValid[blockId]) {
            blockIds.push(blockId);
            readBlockFiles[blockId].open(blockFileName, std::ios::in | std::ios::binary);
        }
    }
    logger->info("Block counts for next step is {}.", blockIds.size());
    for (uint32_t i = 0; i < nThreads; ++i) {
        threads.emplace_back(std::thread(&MSTMerger::buildPairedColorIdEdgesInParallel, this, i,
                                         std::ref(*cqf), std::ref(blockIds), std::ref(readBlockFiles),
                                         std::ref(blockCnt),
                                         std::ref(maxId)));
    }
    for (auto &t : threads) { t.join(); }

    cqf->free();

    logger->info("Total number of kmers observed: {}", numOfKmers);
    num_colorClasses = maxId + 1;
    logger->info("Put edges in each bucket in a sorted list.");
    for (uint32_t i = 0; i < nThreads; ++i) {
        std::string filename = prefix + "/tmp" + std::to_string(i);
        std::ifstream tmp;
        tmp.open(filename, std::ios::in | std::ios::binary);
        uint64_t cnt;
        tmp.read(reinterpret_cast<char *>(&cnt), sizeof(cnt));
        logger->info("file {} has {} edges.", i, cnt);
        std::vector<Edge> edgeList;
        edgeList.resize(cnt);
        tmp.read(reinterpret_cast<char *>(edgeList.data()), sizeof(Edge) * cnt);
        tmp.close();
        std::remove(filename.c_str());
        std::sort(edgeList.begin(), edgeList.end(),
                  [](Edge &e1, Edge &e2) {
                      return e1.n1 == e2.n1 ? e1.n2 < e2.n2 : e1.n1 < e2.n1;
                  });
        edgeList.erase(std::unique(edgeList.begin(), edgeList.end(),
                                   [](Edge &e1, Edge &e2) {
                                       return e1.n1 == e2.n1 and e1.n2 == e2.n2;
                                   }), edgeList.end());
        for (auto &edge: edgeList) {
            edgeBucketList[getBucketId(edge.n1, edge.n2)].push_back(edge);
        }
    }
    for (auto &bucket: edgeBucketList) {
        std::cerr << "before uniqifying: " << bucket.size() << " ";
        std::sort(bucket.begin(), bucket.end(),
                  [](Edge &e1, Edge &e2) {
                      return e1.n1 == e2.n1 ? e1.n2 < e2.n2 : e1.n1 < e2.n1;
                  });
        bucket.erase(std::unique(bucket.begin(), bucket.end(),
                                 [](Edge &e1, Edge &e2) {
                                     return e1.n1 == e2.n1 and e1.n2 == e2.n2;
                                 }), bucket.end());
        std::cerr << "after: " << bucket.size() << "\n";
    }
    logger->info("Done sorting the edges.");

    // Add an edge between each color class ID and node zero
    logger->info("Adding edges from dummy node zero to each color class Id for {} color classes",
                 num_colorClasses);
    zero = static_cast<colorIdType>(num_colorClasses);
    for (colorIdType colorId = 0; colorId < num_colorClasses; colorId++) {
        edgeBucketList[getBucketId(colorId, zero)].push_back(Edge(colorId, zero));
    }
    num_colorClasses++; // zero is now a dummy color class with ID equal to actual num of color classes
    return true;
}

void MSTMerger::writePotentialColorIdEdgesInParallel(uint32_t threadId,
                                                     CQF<KeyObject> &cqf,
                                                     std::queue<uint64_t> &blockIds,
                                                     std::vector<std::ofstream> &blockFiles,
                                                     std::vector<uint64_t> &blockCnt) {
    logger->info("writePotentialColorIdEdgesInParallel {}", threadId);
    uint64_t blockCntr{0};
    std::vector<std::vector<std::pair<colorIdType, uint32_t>>> localBlocks(NUMBlocks);
    for (auto &b : localBlocks) {
        b.reserve(MAX_ALLOWED_BLOCK_SIZE);
    }
    while (true) {
        uint64_t blockId = 0;
        colorMutex.lock();
        if (blockIds.empty()) {
            std::cerr << "\r";
            logger->info("{}: done. Exiting", threadId);
            colorMutex.unlock();
            return;
        }
        blockId = blockIds.front();
        blockIds.pop();
        colorMutex.unlock();
        blockCntr++;
        if (blockId % 10000 == 0)
            std::cerr << "\r" << threadId << "->" << blockId;
        auto nextBlockId = blockId + 1;
        blockId <<= SHIFTBITS;
        nextBlockId <<= SHIFTBITS;
        __uint128_t startPoint = cqf.range() < blockId ? cqf.range() : blockId;
        __uint128_t endPoint = cqf.range() < nextBlockId ? cqf.range() : nextBlockId;
        auto it = cqf.setIteratorLimits(startPoint, endPoint);
        while (!it.reachedHashLimit()) {
            KeyObject keyObject = *it;
            uint64_t curEqId = keyObject.count - 1;
            dna::canonical_kmer n(static_cast<int>(k), keyObject.key);
            // Add an edge between the color class and each of its neighbors' colors in dbg
            for (const auto b : dna::bases) {
                dna::canonical_kmer newNode = n << b;
                if (newNode.val > n.val) {
                    uint64_t idx = (newNode.val & LOWBIT_MASK) >> SHIFTBITS;
                    localBlocks[idx].emplace_back(curEqId, newNode.val);
                    if (localBlocks[idx].size() == MAX_ALLOWED_BLOCK_SIZE) {
                        blockMutex[idx]->lock();
                        blockFiles[idx].write(reinterpret_cast<char *>(localBlocks[idx].data()),
                                              (sizeof(colorIdType) + sizeof(uint32_t)) * MAX_ALLOWED_BLOCK_SIZE);
                        blockCnt[idx] += localBlocks[idx].size();
                        blockMutex[idx]->unlock();
                        localBlocks[idx].clear();
                    }
                }

                newNode = b >> n;
                if (newNode.val > n.val) {
                    uint64_t idx = (newNode.val & LOWBIT_MASK) >> SHIFTBITS;
                    localBlocks[idx].emplace_back(curEqId, newNode.val);
                    if (localBlocks[idx].size() == MAX_ALLOWED_BLOCK_SIZE) {
                        blockMutex[idx]->lock();
                        blockFiles[idx].write(reinterpret_cast<char *>(localBlocks[idx].data()),
                                              (sizeof(colorIdType) + sizeof(uint32_t)) * MAX_ALLOWED_BLOCK_SIZE);
                        blockCnt[idx] += localBlocks[idx].size();
                        blockMutex[idx]->unlock();
                        localBlocks[idx].clear();
                    }
                }
            }
            ++it;
        }
//        if (blockCntr % 10 == 0) {
            for (auto idx = 0; idx < localBlocks.size(); idx++) {
                blockMutex[idx]->lock();
                if (!localBlocks[idx].empty())
                    blockFiles[idx].write(reinterpret_cast<char *>(localBlocks[idx].data()),
                                          (sizeof(colorIdType) + sizeof(uint32_t)) * localBlocks[idx].size());
//                std::cerr << localBlocks[idx].size() << "\n";
                blockCnt[idx] += localBlocks[idx].size();
                blockMutex[idx]->unlock();
                localBlocks[idx].clear();
            }
//        }
    }
}

void MSTMerger::buildPairedColorIdEdgesInParallel(uint32_t threadId,
                                                  CQF<KeyObject> &cqf,
                                                  std::queue<uint64_t> &blockIds,
                                                  std::vector<std::ifstream> &blockFiles,
                                                  std::vector<uint64_t> &blockCnt,
                                                  uint64_t &maxId) {

    logger->info("buildPairedColorIdEdgesInParallel {}", threadId);
    uint64_t blockId = 0;
    uint64_t localMaxId{0}, edgeCntr{0};
    colorMutex.lock();
    if (blockIds.empty()) {
        colorMutex.unlock();
        return;
    }
    blockId = blockIds.front();
    blockIds.pop();
    colorMutex.unlock();
    std::string filename(prefix + "/tmp" + std::to_string(threadId));
    std::ofstream tmpfile(filename, std::ios::out | std::ios::binary);
    tmpfile.write(reinterpret_cast<const char *>(&edgeCntr), sizeof(edgeCntr));
    std::vector<std::pair<colorIdType, uint32_t>> blockBuffer;
    while (true) {
        uint64_t blockKmerCnt{blockCnt[blockId]};
        blockBuffer.clear();
//        logger->info("{}:{}", blockId, blockKmerCnt);
        blockBuffer.resize(blockKmerCnt);
        blockFiles[blockId].read(reinterpret_cast<char *>(blockBuffer.data()),
                                 blockKmerCnt * (sizeof(colorIdType) + sizeof(uint32_t)));
        blockFiles[blockId].close();

        std::string blockFileName = prefix + "/tmpblock" + std::to_string(blockId);
        std::remove(blockFileName.c_str());

        auto nextBlockId = blockId + 1;

        blockId <<= SHIFTBITS;
        nextBlockId <<= SHIFTBITS;
        __uint128_t startPoint = cqf.range() < blockId? cqf.range() : blockId;
        __uint128_t endPoint = cqf.range() < nextBlockId? cqf.range() : nextBlockId;
        std::sort(blockBuffer.begin(), blockBuffer.end(),
                  [](auto &e1, auto &e2) {
                      return e1.second == e2.second ? e1.first < e2.first : e1.second < e2.second;
                  });
        blockBuffer.erase(std::unique(blockBuffer.begin(), blockBuffer.end(),
                                      [](auto &e1, auto &e2) {
                                          return e1.first == e2.first and e1.second == e2.second;
                                      }), blockBuffer.end());
        uint64_t bCntr{0};
        auto tmpEdgeListSize = MAX_ALLOWED_TMP_EDGES / nThreads;
        std::vector<Edge> edgeList;
        edgeList.reserve(tmpEdgeListSize);
        auto it = cqf.setIteratorLimits(startPoint, endPoint);
        KeyObject keyObject = *it;
        dna::canonical_kmer n(static_cast<int>(k), keyObject.key);
        uint32_t curKey{0}, blockKey{0}, blockColorId{0};
        uint64_t curColorId{0};
        bool fetchCQFVal{true}, fetchBlockVal{true};
        while (!it.reachedHashLimit() and bCntr < blockBuffer.size()) {
            if (fetchCQFVal) {
                keyObject = *it;
                n = dna::canonical_kmer(static_cast<int>(k), keyObject.key);
                curKey = static_cast<uint32_t > (n.val & HIGHBIT_MASK);
                curColorId = keyObject.count - 1;
                fetchCQFVal = false;
            }
            if (fetchBlockVal) {
                blockColorId = blockBuffer[bCntr].first;
                blockKey = blockBuffer[bCntr].second;
                fetchBlockVal = false;
            }
            if (curKey < blockKey) {
                //curKey++
                ++it;
                fetchCQFVal = true;
            } else if (blockKey < curKey) {
                //blockKey++
                bCntr++;
                fetchBlockVal = true;
            } else { // curKey == blockKey
                //create and add the new edge
                if (blockColorId != curColorId) {
                    edgeList.emplace_back(blockColorId, curColorId);
                }
                //curKey++
                ++it;
                fetchCQFVal = true;
                //blockKey++
                bCntr++;
                fetchBlockVal = true;
                localMaxId = curColorId > localMaxId ? (curColorId > blockColorId ? curColorId : blockColorId) :
                             (blockColorId > localMaxId ? blockColorId : localMaxId);
            }
            if (edgeList.size() >= tmpEdgeListSize/* and colorMutex.try_lock()*/) {
                tmpfile.write(reinterpret_cast<const char *>(edgeList.data()), sizeof(Edge) * edgeList.size());
                edgeCntr += edgeList.size();
                edgeList.clear();
            }
        }
        tmpfile.write(reinterpret_cast<const char *>(edgeList.data()), sizeof(Edge) * edgeList.size());
        edgeCntr += edgeList.size();
        colorMutex.lock();
        if (blockIds.empty()) {
            maxId = localMaxId > maxId ? localMaxId : maxId;
            colorMutex.unlock();
            std::cerr << "\r";
            logger->info("Thread {}: added {} edges", threadId, edgeCntr);
            tmpfile.seekp(0);
            tmpfile.write(reinterpret_cast<const char *>(&edgeCntr), sizeof(edgeCntr));
            tmpfile.close();
            return;
        }
        blockId = blockIds.front();
        blockIds.pop();
        colorMutex.unlock();
    }
}


/**
 * loads the color class table in parts
 * calculate the hamming distance between the color bitvectors fetched from color class table
 * for each pair of color IDs
 * having w buckets where w is the maximum possible weight (number of experiments)
 * put the pair in its corresponding bucket based on the hamming distance value (weight)
 * @return true if successful
 */
bool MSTMerger::calculateMSTBasedWeights() {

    mst1 = new MSTQuery(prefix1, k, k, numOfFirstMantisSamples, logger);
    mst2 = new MSTQuery(prefix2, k, k, secondMantisSamples, logger);

    nonstd::optional<uint64_t> dummy{nonstd::nullopt};

    QueryStats dummyStats1, dummyStats2;

    logger->info("loaded the two msts with k={}. MST sizes are {}, {} respectively.", k, mst1->parentbv.size(),
                 mst2->parentbv.size());
    std::ifstream cp(prefix + "newID2oldIDs");
    uint64_t cnt, cIdx, c1, c2;
    cp.read(reinterpret_cast<char *>(&cnt), sizeof(cnt));
    logger->info("# of color classes based on count of colorPairs: {}", cnt);
    colorPairs.resize(cnt);
    for (auto i = 0; i < cnt; i++) {
        cp.read(reinterpret_cast<char *>(&cIdx), sizeof(cIdx));
        cp.read(reinterpret_cast<char *>(&c1), sizeof(c1));
        cp.read(reinterpret_cast<char *>(&c2), sizeof(c2));
        c1 = c1 == 0 ? mst1->parentbv.size() - 1 : c1 - 1;
        c2 = c2 == 0 ? mst2->parentbv.size() - 1 : c2 - 1;
        //std::cerr << cIdx << " " << n1s << " " << n2s << "\n";
        colorPairs[cIdx] = std::make_pair(c1, c2);
    }
    cp.close();

    logger->info("Splitting the edges into edges for MST1 and MST2 for {} eqclass buckets.", eqclass_files.size());
    uint64_t numEdges = 0;
    weightBuckets.resize(numSamples);
    for (auto &qf : queryStats1) {
        qf.numSamples = numSamples;
        qf.trySample = true;
    }
    for (auto &qf : queryStats2) {
        qf.numSamples = numSamples;
        qf.trySample = true;
    }
    colorIdType mst1Zero = mst1->parentbv.size() - 1, mst2Zero = mst2->parentbv.size() - 1;
    std::vector<std::pair<colorIdType, weightType>> edge1list;
    std::vector<std::pair<colorIdType, weightType>> edge2list;
    std::vector<uint32_t> srcStartIdx1(mst1->parentbv.size(), 0), srcStartIdx2(mst2->parentbv.size(), 0);


    // total merged edges can be a good estimate to reserve a vector per mst
    uint64_t totalEdgeCnt{0};
    for (auto i = 0; i < eqclass_files.size(); i++) {
        for (auto j = i; j < eqclass_files.size(); j++) {
            auto &edgeBucket = edgeBucketList[i * num_of_ccBuffers + j];
            totalEdgeCnt += edgeBucket.size();
        }
    }
    logger->info("Total # of edges: {}", totalEdgeCnt);

    // lambda function to get sorted list of unique edges for each mantis.
    auto edgeSplitter = [=](std::vector<uint32_t> &srcStartIdx,
                            std::vector<std::pair<colorIdType, weightType>> &destWeightList,
                            colorIdType mstZero,
                            bool isFirst) {
        std::vector<Edge> edgeList;
        edgeList.reserve(totalEdgeCnt);
        // put all of the edges for one of the merged mantises in a vector
        for (auto i = 0; i < eqclass_files.size(); i++) {
            for (auto j = i; j < eqclass_files.size(); j++) {
                auto &edgeBucket = edgeBucketList[i * num_of_ccBuffers + j];
                for (auto &edge : edgeBucket) {
                    auto n1 = edge.n1 == zero ? mstZero : (isFirst ? colorPairs[edge.n1].first
                                                                   : colorPairs[edge.n1].second);
                    auto n2 = edge.n2 == zero ? mstZero : (isFirst ? colorPairs[edge.n2].first
                                                                   : colorPairs[edge.n2].second);
                    if (n1 != n2)
                        edgeList.emplace_back(n1, n2);
                }
                std::cerr << "\rmantis" << (isFirst ? "1" : "2") << " edgelist for " << i << "," << j << "  ";
            }
        }
        std::cerr << "\r";
        // sort - unique the vector
        std::sort(edgeList.begin(), edgeList.end(),
                  [](Edge &e1, Edge &e2) {
                      return e1.n1 == e2.n1 ? e1.n2 < e2.n2 : e1.n1 < e2.n1;
                  });
        edgeList.erase(std::unique(edgeList.begin(), edgeList.end(),
                                   [](Edge &e1, Edge &e2) {
                                       return e1.n1 == e2.n1 and e1.n2 == e2.n2;
                                   }), edgeList.end());

        // move all the unique edges for one mantis to the designed low-memory data structure
        destWeightList.resize(edgeList.size());
        uint64_t destIdx{0}, prevCnt{0};
        for (auto &e : edgeList) {
            srcStartIdx[e.n1]++;
            destWeightList[destIdx++] = std::make_pair(e.n2, 0);
        }
        // accumulate the out-degrees per sorted source edges to get to source start index
        // now each value in srcStartIdx is pointing to the start of the next source in the edge-weight vec.
        for (auto &outCnt : srcStartIdx) {
            outCnt += prevCnt;
            prevCnt = outCnt;
        }
    };
    // call the lambda function per mantis
    edgeSplitter(srcStartIdx1, edge1list, mst1Zero, true);
    edgeSplitter(srcStartIdx2, edge2list, mst2Zero, false);
/*
    std::cerr << "\n";
    for (auto i = 0; i < 10; i++) {
        std::cerr << i << ":" << srcStartIdx1[i] << "\n";
        auto j = i == 0? 0 : srcStartIdx1[i-1];
        while (j < srcStartIdx1[i]) {
            std::cerr << j << ":" << edge1list[j].first << " ";
            j++;
        }
        std::cerr << "\n";
    }*/
    logger->info("num of edges in first mantis: {}", edge1list.size());
    logger->info("num of edges in second mantis: {}", edge2list.size());

    std::vector<colorIdType> colorsInCache;
    planCaching(mst1, edge1list, srcStartIdx1, colorsInCache);
    logger->info("fixed cache size for mst1 is : {}", colorsInCache.size());
    // fillout fixed_cache1
    // walking in reverse order on colorsInCache is intentional
    // because as we've filled out the colorsInCache vector in a post-order traversal of MST,
    // if we walk in from end to the beginning we can always guarantee that we've already put the ancestors
    // for the new colors we want to put in the fixed_cache and that will make the whole color bv construction
    // FASTER!!
    for (int64_t idx = colorsInCache.size() - 1; idx >= 0; idx--) {
        auto setbits = mst1->buildColor(colorsInCache[idx], dummyStats1, &lru_cache1[0], nullptr, &fixed_cache1, dummy);
        fixed_cache1[colorsInCache[idx]] = setbits;
    }
    colorsInCache.clear();

    logger->info("Done filling the fixed cache for mst1. Calling multi-threaded MSTBasedHammingDist .. ");
    std::vector<std::thread> threads;
    for (uint32_t t = 0; t < nThreads; ++t) {
        threads.emplace_back(std::thread(&MSTMerger::calcMSTHammingDistInParallel, this, t,
                                         std::ref(edge1list),
                                         std::ref(srcStartIdx1),
                                         mst1,
                                         std::ref(lru_cache1),
                                         std::ref(queryStats1),
                                         std::ref(fixed_cache1),
                                         numOfFirstMantisSamples));
    }
    for (auto &t : threads) { t.join(); }
    threads.clear();
    logger->info("Done calculating Dist for mst1");


    planCaching(mst2, edge2list, srcStartIdx2, colorsInCache);
    logger->info("fixed cache size for mst2 is : {}", colorsInCache.size());
//    std::exit(3);
    //fillout fixed_cache2
    for (int64_t idx = colorsInCache.size() - 1; idx >= 0; idx--) {
        auto setbits = mst2->buildColor(colorsInCache[idx], dummyStats2, &lru_cache2[0], nullptr, &fixed_cache2, dummy);
        fixed_cache2[colorsInCache[idx]] = setbits;
    }

    logger->info("Done filling the fixed cache for mst2. Calling multi-threaded MSTBasedHammingDist .. ");
    for (uint32_t t = 0; t < nThreads; ++t) {
        threads.emplace_back(std::thread(&MSTMerger::calcMSTHammingDistInParallel, this, t,
                                         std::ref(edge2list),
                                         std::ref(srcStartIdx2),
                                         mst2,
                                         std::ref(lru_cache2),
                                         std::ref(queryStats2),
                                         std::ref(fixed_cache2),
                                         secondMantisSamples));
    }
    for (auto &t : threads) { t.join(); }
    colorsInCache.clear();
    logger->info("Done calculating Dist for mst2");
    uint64_t cntr{0};

    //abstract out the part repeated for mst1 and mst2 as a lambda function
    auto findWeight = [=](Edge &edge, std::vector<std::pair<colorIdType, weightType>> &edgeList,
                          std::vector<uint32_t> &srcStartIdx, colorIdType mstZero, bool isFirst) {
//            std::cerr << edge.n1 << " " << edge.n2 << "\n";
//            std::cerr << mstZero << "\n";
        auto n1 = edge.n1 == zero ? mstZero : (isFirst ? colorPairs[edge.n1].first : colorPairs[edge.n1].second);
        auto n2 = edge.n2 == zero ? mstZero : (isFirst ? colorPairs[edge.n2].first : colorPairs[edge.n2].second);
//            std::cerr << n1 << " " << n2 << " ";
        if (n1 == n2) return static_cast<weightType>(0);
        if (n1 > n2) std::swap(n1, n2);
        weightType w{0};
        auto srcStart = n1 == 0 ? 0 : srcStartIdx[n1 - 1];
        if (srcStart < edgeList.size()) {
            auto srcEnd = srcStartIdx[n1];
            if (n2 == mstZero) { // zero is the biggest color ID, and hence the last in a sorted list
                if ((*(edgeList.begin() + srcEnd - 1)).first != n2) {
                    std::cerr << "!!NOOOOOO!\n";
                    std::cerr << n1 << " " << srcStart << " " << srcEnd << " " << n2 << " "
                              << (*(edgeList.begin() + srcEnd - 1)).first << "\n";
                    std::exit(3);
                }
                w = (*(edgeList.begin() + srcEnd - 1)).second; // look at the last elm. for n1
            } else {
                auto wItr = std::lower_bound(edgeList.begin() + srcStart, edgeList.begin() + srcEnd,
                                             std::make_pair(n2, 0),
                                             [](auto &v1, auto &v2) {
                                                 return v1.first < v2.first;
                                             });
                if (wItr == edgeList.begin() + srcEnd) {
                    std::cerr << "Couldn't find the element in the vector. Should not happen.\n";
                    std::exit(3);
                }
                if ((*wItr).first != n2) {
                    std::cerr << "NOOOOOO!\n";
                    std::cerr << n1 << " " << srcStart << " " << srcEnd << " " << n2 << " " << (*wItr).first << "\n";
                    std::exit(3);
                }
                w = (*wItr).second;
            }
        }
//        std::cerr << n1 << " " << n2 << " " << w << "\n";
        return w;
    };
    for (auto i = 0; i < eqclass_files.size(); i++) {
        for (auto j = i; j < eqclass_files.size(); j++) {
            auto &edgeBucket = edgeBucketList[i * num_of_ccBuffers + j];
            for (auto &edge : edgeBucket) {
//                std::cerr << edge.n1 << " " << edge.n2 << "\n";
                auto w1 = findWeight(edge, edge1list, srcStartIdx1, mst1Zero, true);
                auto w2 = findWeight(edge, edge2list, srcStartIdx2, mst2Zero, false);
//                std::cerr << w1 << " " << w2 << "\n";
                weightBuckets[w1 + w2 - 1].push_back(edge);
                if (++cntr % 10000000 == 0)
                    std::cerr << "\r" << cntr << " out of " << totalEdgeCnt;
            }
            edgeBucket.clear();
        }
    }
    std::cerr << "\r";
    edgeBucketList.clear();
    srcStartIdx1.clear();
    srcStartIdx2.clear();
    edge1list.clear();
    edge2list.clear();
    std::cerr << "\r";
    logger->info("Calculated the weight for the edges");
    /*logger->info("Writing the distributions down ..");
    std::ofstream cacheHeight(prefix + "/cacheHeight.dist");
    std::ofstream cacheWeight(prefix + "/cacheWeight.dist");
    std::ofstream noCacheHeight(prefix + "/noCacheHeight.dist");
    std::ofstream noCacheWeight(prefix + "/noCacheWeight.dist");
    for (auto &qs : queryStats1) {
        for (auto v : qs.heightDist) cacheHeight << "1 " << v << "\n";
        for (auto v : qs.weightDist) cacheWeight << "1 " << v << "\n";
        for (auto v : qs.noCache_heightDist) noCacheHeight << "1 " << v << "\n";
        for (auto v : qs.noCache_weightDist) noCacheWeight << "1 " << v << "\n";
    }
    for (auto &qs : queryStats2) {
        for (auto v : qs.heightDist) cacheHeight << "2 " << v << "\n";
        for (auto v : qs.weightDist) cacheWeight << "2 " << v << "\n";
        for (auto v : qs.noCache_heightDist) noCacheHeight << "2 " << v << "\n";
        for (auto v : qs.noCache_weightDist) noCacheWeight << "2 " << v << "\n";
    }
    cacheHeight.close();
    cacheWeight.close();
    noCacheHeight.close();
    noCacheWeight.close();
    logger->info("Done Writing the distributions down");*/
    return true;
}


void MSTMerger::calcMSTHammingDistInParallel(uint32_t i,
                                             std::vector<std::pair<colorIdType, weightType>> &edgeList,
                                             std::vector<uint32_t> &srcStarts,
                                             MSTQuery *mst,
                                             std::vector<LRUCacheMap> &lru_cache,
                                             std::vector<QueryStats> &queryStats,
                                             std::unordered_map<uint64_t, std::vector<uint64_t>> &fixed_cache,
                                             uint32_t numSamples) {
    std::vector<uint64_t> srcBV;
    uint64_t s = 0, e;
    // If the list contains less than a hundred edges, don't bother with multi-threading and
    // just run the first thread
    if (edgeList.size() < 100/* or isMSTBased*/) {
        if (i == 0) {
            e = edgeList.size();
        } else {
            e = 0;
        }
    } else {
        s = edgeList.size() * i / nThreads;
        e = edgeList.size() * (i + 1) / nThreads;
    }
    auto start = std::upper_bound(srcStarts.begin(), srcStarts.end(), s);
    colorIdType n1 = static_cast<colorIdType>(start - srcStarts.begin());
    auto itrStart = edgeList.begin() + s;
    auto itrEnd = edgeList.begin() + e;
    uint64_t cntr{0}, edgeCntr{s};
    std::cerr << "\r";
//    logger->info("initiated thread {} from {} to {} with {} edges", i, s, e, e-s);
    for (auto edge = itrStart; edge != itrEnd; edge++) {
        if (srcStarts[n1] == edgeCntr) {
            n1++;
        }
        colorIdType n2 = edge->first;
        if (n1 != n2) {
            auto w = mstBasedHammingDist(n1, n2,
                                         mst, lru_cache[i], srcBV, queryStats[i],
                                         fixed_cache);
            edge->second = static_cast<weightType>(w);
        }
        if (++cntr % 1000000 == 0) {
            std::cerr << "\rthread " << i << ": " << cntr << " edges out of " << e - s << "   ";
        }
        edgeCntr++;
    }
    std::cerr << "\r";
}

/**
 * Finds Minimum Spanning Forest of color graph using Kruskal Algorithm
 *
 * The algorithm's basic implementation taken from
 * https://www.geeksforgeeks.org/kruskals-minimum-spanning-tree-using-stl-in-c/
 * @return List of connected components in the Minimum Spanning Forest
 */
DisjointSets MSTMerger::kruskalMSF() {
    uint32_t bucketCnt = numSamples;
    mst.resize(num_colorClasses);
    // Create disjoint sets
    DisjointSets ds(num_colorClasses);

    uint64_t edgeCntr{0}, selectedEdgeCntr{0};
    uint32_t w{0};

    // Iterate through all sorted edges
    for (uint32_t bucketCntr = 0; bucketCntr < bucketCnt; bucketCntr++) {
        uint32_t edgeIdxInBucket = 0;
        w = bucketCntr + 1;
        for (auto &it : weightBuckets[bucketCntr]) {
            colorIdType u = it.n1;
            colorIdType v = it.n2;
            colorIdType root_of_u = ds.find(u);
            colorIdType root_of_v = ds.find(v);

            // Check if the selected edge is causing a cycle or not
            // (A cycle is induced if u and v belong to the same set)
            if (root_of_u != root_of_v) {
                // Merge two sets
                ds.merge(root_of_u, root_of_v, w);
                // Current edge will be in the MST
                mst[u].emplace_back(v, w);
                mst[v].emplace_back(u, w);
                mstTotalWeight += w;
                selectedEdgeCntr++;
            }
            edgeCntr++;
            if (edgeCntr % 1000000 == 0) {
                std::cerr << "\r" << edgeCntr << " edges processed and "
                          << selectedEdgeCntr << " were selected";
            }
            edgeIdxInBucket++;
        }
        weightBuckets[bucketCntr].clear();
    }
    std::cerr << "\r";
    mstTotalWeight++;//1 empty slot for root (zero)
    logger->info("MST Construction finished:"
                 "\n\t# of graph edges: {}"
                 "\n\t# of merges (mst edges): {}"
                 "\n\tmst weight sum: {}",
                 edgeCntr, selectedEdgeCntr, mstTotalWeight);
    return ds;
}

/**
 * calls kruskal algorithm to build an MST of the color graph
 * goes over the MST and fills in the int-vectors parentbv, bbv, and deltabv
 * serializes these three int-vectors as the encoding of color classes
 * @return true if encoding and serializing the DS is successful
 */
bool MSTMerger::encodeColorClassUsingMST() {
    // build mst of color class graph
    kruskalMSF();

    uint64_t nodeCntr{0};
    // encode the color classes using mst
    logger->info("Filling ParentBV...");
    sdsl::int_vector<> parentbv(num_colorClasses, 0, ceil(log2(num_colorClasses)));
    // create and fill the deltabv and boundarybv data structures
    sdsl::bit_vector bbv(mstTotalWeight, 0);
    {// putting weightbv inside the scope so its memory is freed after we're done with it
        sdsl::int_vector<> weightbv(num_colorClasses, 0, ceil(log2(numSamples)));
        sdsl::bit_vector visited(num_colorClasses, 0);
        bool check = false;
        std::queue<colorIdType> q;
        q.push(zero); // Root of the tree is zero
        parentbv[zero] = zero; // and it's its own parent (has no parent)
        weightbv[zero] = 1; // adding a dummy weight for a dummy node
        while (!q.empty()) {
            colorIdType parent = q.front();
            q.pop();
            for (auto &neighbor :mst[parent]) {
                if (!visited[neighbor.first]) {
                    parentbv[neighbor.first] = parent;
                    weightbv[neighbor.first] = neighbor.second;
                    q.push(neighbor.first);
                }
            }
            visited[parent] = 1;
            nodeCntr++; // just a counter for the log
            if (nodeCntr % 10000000 == 0) {
                std::cerr << "\rset parent of " << nodeCntr << " ccs";
            }
        }

        std::cerr << "\r";
        // filling bbv
        // resize bbv
        logger->info("Filling BBV...");
        uint64_t deltaOffset{0};
        for (uint64_t i = 0; i < num_colorClasses; i++) {
            deltaOffset += static_cast<uint64_t>(weightbv[i]);
            bbv[deltaOffset - 1] = 1;
        }
    }
    std::cerr << "\r";
    // fill in deltabv
    logger->info("Filling DeltaBV...");
    sdsl::int_vector<> deltabv(mstTotalWeight, 0, ceil(log2(numSamples)));
    sdsl::bit_vector::select_1_type sbbv = sdsl::bit_vector::select_1_type(&bbv);
    std::vector<std::thread> threads;
    for (uint32_t t = 0; t < nThreads; ++t) {
        threads.emplace_back(std::thread(&MSTMerger::calcDeltasInParallel, this,
                                         t, 0, 0,
                                         std::ref(parentbv), std::ref(deltabv), std::ref(sbbv), true));
    }
    for (auto &t : threads) { t.join(); }
    std::cerr << "\r";

    logger->info("Serializing data structures parentbv, deltabv, & bbv...");
    sdsl::store_to_file(parentbv, std::string(prefix + mantis::PARENTBV_FILE));
    sdsl::store_to_file(deltabv, std::string(prefix + mantis::DELTABV_FILE));
    sdsl::store_to_file(bbv, std::string(prefix + mantis::BOUNDARYBV_FILE));
    logger->info("Done Serializing.");
    return true;
}

void MSTMerger::calcDeltasInParallel(uint32_t threadID, uint64_t cbvID1, uint64_t cbvID2,
                                     sdsl::int_vector<> &parentbv, sdsl::int_vector<> &deltabv,
                                     sdsl::bit_vector::select_1_type &sbbv,
                                     bool isMSTBased) {

    struct Delta {
        uint64_t startingOffset{0};
        std::vector<uint32_t> deltaVals;

        Delta() = default;

        Delta(uint64_t so) {
            startingOffset = so;
        }
    };
    std::vector<Delta> deltas;

    colorIdType s = parentbv.size() * threadID / nThreads;
    colorIdType e = parentbv.size() * (threadID + 1) / nThreads;

    colorIdType mst1Zero = mst1->parentbv.size() - 1, mst2Zero = mst2->parentbv.size() - 1;
    for (colorIdType p = s; p < e; p++) {
        auto deltaOffset = (p > 0) ? (sbbv(p) + 1) : 0;
        deltas.push_back(deltaOffset);
        auto n1s = p == zero ? std::make_pair(mst1Zero, mst2Zero) : colorPairs[p];
        auto n2s = parentbv[p] == zero ? std::make_pair(mst1Zero, mst2Zero) : colorPairs[parentbv[p]];
        auto firstDelta = getMSTBasedDeltaList(n1s.first, n2s.first, lru_cache1[threadID], true,
                                               queryStats1[threadID]);
        auto secondDelta = getMSTBasedDeltaList(n1s.second, n2s.second, lru_cache2[threadID], false,
                                                queryStats2[threadID]);
        deltas.back().deltaVals = firstDelta;
        for (auto &v : secondDelta) {
            v += numOfFirstMantisSamples;
            deltas.back().deltaVals.push_back(v);
        }
    }

    colorMutex.lock();
    for (auto &v : deltas) {
        for (auto cntr = 0; cntr < v.deltaVals.size(); cntr++)
            deltabv[v.startingOffset + cntr] = v.deltaVals[cntr];
    }
    colorMutex.unlock();

}

/**
 * finds the neighbors of each kmer in the cqf,
 * and adds an edge of the element's colorId and its neighbor's
 * @param cqf (required to query for existence of neighbors)
 * @param it iterator to the elements of cqf
 */
void MSTMerger::findNeighborEdges(CQF<KeyObject> &cqf, KeyObject &keyobj, std::vector<Edge> &edgeList) {
    dna::canonical_kmer curr_node(static_cast<int>(k), keyobj.key);
    workItem cur = {curr_node, static_cast<colorIdType>(keyobj.count - 1)};
    uint64_t neighborCnt{0};
    for (auto &nei : neighbors(cqf, cur)) {
        neighborCnt++;
        if (cur.colorId < nei.colorId) {
            Edge e(static_cast<colorIdType>(cur.colorId), static_cast<colorIdType>(nei.colorId));
            edgeList.push_back(e);
        }
    }
}

/**
 * finds neighbors of a node in cqf
 * @param cqf
 * @param n : work_item containing node and colorId (colorId will be filled)
 * @return set of neighbors for current node n and their colorIds
 */
std::set<workItem> MSTMerger::neighbors(CQF<KeyObject> &cqf, workItem n) {
    std::set<workItem> result;
    for (const auto b : dna::bases) {
        uint64_t eqid = 0;
        if (exists(cqf, n.node << b, eqid)) {
            if (eqid != n.colorId)
                result.insert(workItem(n.node << b, eqid));
        }
        if (exists(cqf, b >> n.node, eqid)) {
            if (eqid != n.colorId)
                result.insert(workItem(b >> n.node, eqid));
        }
    }
    return result;
}

/**
 * searches for a kmer in cqf and returns the correct colorId if found
 * which is cqf count value - 1
 * @param cqf
 * @param e : search canonical kmer
 * @param eqid : reference to eqid that'll be set
 * @return true if eqid is found
 */
bool MSTMerger::exists(CQF<KeyObject> &cqf, dna::canonical_kmer e, uint64_t &eqid) {
    KeyObject key(e.val, 0, 0);
    auto eqidtmp = cqf.query(key, QF_NO_LOCK /*QF_KEY_IS_HASH | QF_NO_LOCK*/);
    if (eqidtmp) {
        eqid = eqidtmp - 1;
        return true;
    }
    return false;
}

void MSTMerger::buildMSTBasedColor(uint64_t eqid,
                                   MSTQuery *mst,
                                   LRUCacheMap &lru_cache,
                                   std::vector<uint64_t> &eq,
                                   QueryStats &queryStats,
                                   std::unordered_map<uint64_t, std::vector<uint64_t>> &fixed_cache) {
    RankScores rs(1);

    nonstd::optional<uint64_t> dummy{nonstd::nullopt};

//    auto eq_ptr = lru_cache.lookup_ts(eqid);
//    if (eq_ptr) {
    if (fixed_cache.find(eqid) != fixed_cache.end()) {
//        std::cerr << "happens! ";
        eq = fixed_cache[eqid];
        queryStats.cacheCntr++;
//        queryStats.heightDist.push_back(0);
//        queryStats.weightDist.push_back(0);
    } else if (lru_cache.contains(eqid)) {
        eq = lru_cache[eqid];//.get(eqclass_id);
//        eq = *eq_ptr;
        queryStats.cacheCntr++;
//        queryStats.heightDist.push_back(0);
//        queryStats.weightDist.push_back(0);
    } else {
        nonstd::optional<uint64_t> toDecode{nonstd::nullopt};
        queryStats.noCacheCntr++;
        queryStats.trySample = (queryStats.noCacheCntr % 20 == 0);
        toDecode.reset();
        eq = mst->buildColor(eqid, queryStats, &lru_cache, &rs, &fixed_cache, toDecode);
//        auto sp = std::make_shared<std::vector<uint64_t>>(eq);
        lru_cache.emplace(eqid, eq);
//        lru_cache.emplace_ts(eqid, sp);
        if (queryStats.trySample and toDecode and fixed_cache.find(*toDecode) == fixed_cache.end()) {
            auto s = mst->buildColor(*toDecode, queryStats, nullptr, nullptr, &fixed_cache, dummy);
            lru_cache.emplace(*toDecode, s);
//            auto sp1 = std::make_shared<std::vector<uint64_t>>(s);
//            lru_cache.emplace_ts(*toDecode, sp1);
        }
    }
}

uint64_t MSTMerger::mstBasedHammingDist(uint64_t eqid1,
                                        uint64_t eqid2,
                                        MSTQuery *mst,
                                        LRUCacheMap &lru_cache,
                                        std::vector<uint64_t> &srcEq,
                                        QueryStats &queryStats,
                                        std::unordered_map<uint64_t, std::vector<uint64_t>> &fixed_cache) {

    uint64_t dist{0};
    std::vector<uint64_t> eq1, eq2;

    buildMSTBasedColor(eqid1, mst, lru_cache, eq1, queryStats, fixed_cache);
    // fetch the second color ID's BV
    buildMSTBasedColor(eqid2, mst, lru_cache, eq2, queryStats, fixed_cache);


    /// calc distance
    auto i{0}, j{0};
    while (i != eq1.size() and j != eq2.size()) {
        if (eq1[i] == eq2[j]) {
            i++;
            j++;
        } else if (eq1[i] < eq2[j]) {
            i++;
            dist++;
        } else {
            j++;
            dist++;
        }
    }
    if (i != eq1.size()) dist += (eq1.size() - i);
    if (j != eq2.size()) dist += (eq2.size() - j);
    return dist;
}

std::vector<uint32_t> MSTMerger::getMSTBasedDeltaList(uint64_t eqid1, uint64_t eqid2, LRUCacheMap &lru_cache,
                                                      bool isFirst, QueryStats &queryStats) {
    std::vector<uint32_t> res;
    if (eqid1 == eqid2) return res;
    std::vector<uint64_t> eq1, eq2;
    if (isFirst) {
        buildMSTBasedColor(eqid1, mst1, lru_cache, eq1, queryStats, fixed_cache1);
        buildMSTBasedColor(eqid2, mst1, lru_cache, eq2, queryStats, fixed_cache1);
    } else {
        buildMSTBasedColor(eqid1, mst2, lru_cache, eq1, queryStats, fixed_cache2);
        buildMSTBasedColor(eqid2, mst2, lru_cache, eq2, queryStats, fixed_cache2);
    }
    /// calc delta
    auto i{0}, j{0};
    while (i != eq1.size() and j != eq2.size()) {
        if (eq1[i] == eq2[j]) {
            i++;
            j++;
        } else if (eq1[i] < eq2[j]) {
            res.push_back(eq1[i]);
            i++;
        } else {
            res.push_back(eq2[j]);
            j++;
        }
    }
    while (i != eq1.size()) {
        res.push_back(eq1[i]);
        i++;
    }
    while (j != eq2.size()) {
        res.push_back(eq2[j]);
        j++;
    }

    return res; // rely on c++ optimization
}


/**
 * calculates the edge corresponding bucket id c1 <= c2
 * @param c1 first colorId
 * @param c2 second colorId
 * @return bucket id
 */
inline uint64_t MSTMerger::getBucketId(uint64_t c1, uint64_t c2) {
    if (c1 == zero or c1 > c2) {
        std::swap(c1, c2);
    }
    uint64_t cb1 = c1 / mantis::NUM_BV_BUFFER;
    uint64_t cb2 = c2 / mantis::NUM_BV_BUFFER;
    if (c2 == zero) // return the corresponding buffer for the non-zero colorId
        return cb1 * num_of_ccBuffers + cb1;
    return cb1 * num_of_ccBuffers + cb2;
}

void MSTMerger::planCaching(MSTQuery *mst,
                            std::vector<std::pair<colorIdType, weightType>> &edges,
                            std::vector<uint32_t> &outDegrees,
                            std::vector<colorIdType> &colorsInCache) {
    std::vector<Cost> mstCost(mst->parentbv.size());

    logger->info("In planner ..");
    // setting local edge costs
    uint64_t src{0}, edgeCntr{0};
    for (auto &edge: edges) {
        if (edgeCntr == outDegrees[src]) {
            src++;
        }
        mstCost[src].numQueries++;
        mstCost[edge.first].numQueries++;
        edgeCntr++;
    }
    logger->info("Done setting the local costs");

    std::vector<std::vector<colorIdType>> children(mst->parentbv.size());
    for (uint64_t i = 0; i < mst->parentbv.size() - 1; i++) {
        children[mst->parentbv[i]].push_back(i);
    }
    logger->info("Done creating the parent->children map");
    // recursive planner
    logger->info("Calling the recursive planner for {} nodes", mst->parentbv.size());
    uint64_t cntr{0};
    planRecursively(mst->parentbv.size() - 1, children, mstCost, colorsInCache, cntr);
    std::cerr << "\r";
}

void MSTMerger::planRecursively(uint64_t nodeId,
                                std::vector<std::vector<colorIdType>> &children,
                                std::vector<Cost> &mstCost,
                                std::vector<colorIdType> &colorsInCache,
                                uint64_t &cntr) {

    uint64_t avgCostThreshold = 16;

//    std::cerr << "\rsize for node " << nodeId << " " << children[nodeId].size();
    for (auto c : children[nodeId]) {
        planRecursively(c, children, mstCost, colorsInCache, cntr);
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
            childWithMaxAvg = 0;
            maxChildAvg = 0;
        }

        numQueries = mstCost[nodeId].numQueries;
        numSteps = mstCost[nodeId].numSteps;
        for (auto c : children[nodeId]) {
            if (localCache.find(c) == localCache.end()) {
                numQueries += mstCost[c].numQueries;
                // steps in parent += steps in children + 1 for each child
                numSteps += mstCost[c].numSteps + mstCost[c].numQueries;
                float currChildAvg =
                        mstCost[c].numQueries == 0 ? 0 : (float) mstCost[c].numSteps / (float) mstCost[c].numQueries;
                if (mstCost[c].numSteps != 0 and
                    (currChildAvg > maxChildAvg or
                     (currChildAvg == maxChildAvg and mstCost[childWithMaxAvg].numQueries < mstCost[c].numQueries))) {
                    maxChildAvg = currChildAvg;
                    childWithMaxAvg = c;
                }
            }
        }
    } while (numQueries != 0 and (float) numSteps / (float) numQueries > avgCostThreshold);
    mstCost[nodeId].numQueries = numQueries;
    mstCost[nodeId].numSteps = numSteps;
    for (auto c : localCache) {
        colorsInCache.push_back(c);
    }
    if (++cntr % 1000000 == 0)
        std::cerr << "\r" << cntr++;
}
