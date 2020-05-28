#include <memory>

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
#include <zconf.h>

#include "MantisFS.h"
#include "mstMerger.h"
#include "ProgOpts.h"


#define BITMASK(nbits) ((nbits) == 64 ? 0xffffffffffffffff : (1ULL << (nbits)) - 1ULL)

static constexpr uint64_t MAX_ALLOWED_TMP_EDGES{31250000}; // Up to 0.5G
static constexpr uint64_t MAX_ALLOWED_TMP_EDGES_IN_FILE{130000000}; // Up to 2G
static constexpr uint64_t SHIFTBITS{32};
static constexpr uint64_t HIGHBIT_MASK{(1ULL << 32) - 1ULL};
static constexpr uint64_t LOWBIT_MASK{std::numeric_limits<uint64_t>::max() - HIGHBIT_MASK};
static constexpr uint64_t CONSTNUMBlocks{1 << 16};

MSTMerger::MSTMerger(std::string prefixIn, spdlog::logger *loggerIn, uint32_t numThreads,
                     std::string prefixIn1, std::string prefixIn2) :
                prefix(std::move(prefixIn)),
                prefix1(std::move(prefixIn1)), prefix2(std::move(prefixIn2)),
                nThreads(numThreads), numBlocks(CONSTNUMBlocks) {
    logger = loggerIn;//.get();

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

    queryStats1.resize(nThreads);
    queryStats2.resize(nThreads);
    for (uint64_t t = 0; t < nThreads; t++) {
        lru_cache1.emplace_back(1000);
        lru_cache2.emplace_back(1000);
    }

    std::string sample_file = prefix1 + mantis::SAMPLEID_FILE;
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

    numCCPerBuffer = mantis::BV_BUF_LEN / numSamples;

    logger->info("# of experiments: {}", numSamples);
    logger->info("# of threads={}, # of ccs per buffer={}, # of edges={}",
            numThreads, numCCPerBuffer, num_edges);
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
    auto t_start = time(nullptr);
    logger->info ("Merging the two MSTs..");

    buildEdgeSets();
    calculateMSTBasedWeights();
    encodeColorClassUsingMST();

//    std::string cmd = "rm " + prefix + "newID2oldIDs";
//    system(cmd.c_str());
    auto t_end = time(nullptr);
    logger->info("MST merge completed in {} s.", t_end - t_start);
}

/**
 * iterates over all elements of CQF,
 * find all the existing neighbors, and build a color graph based on that
 * @return pair of maximum observed colorId and total observed edges (including duplicates)
 */
bool MSTMerger::buildEdgeSets() {

    std::vector<uint64_t> cnts(nThreads, 0);
    uint64_t maxId{0}, numOfKmers{0};


    std::vector<std::pair<uint64_t, uint64_t>> tmpEdges;
    tmpEdges.reserve(MAX_ALLOWED_TMP_EDGES_IN_FILE);
    std::vector<std::string> cqfBlocks = mantis::fs::GetFilesExt(prefix.c_str(), mantis::CQF_FILE);

    for (uint64_t c = 0; c < cqfBlocks.size(); c++) {
        logger->info("Reading colored dbg from disk...");
//        std::cerr << cqfBlocks[c] << "\n";
        std::string cqf_file(cqfBlocks[c]);
        CQF<KeyObject> cqf(cqf_file, CQF_FREAD);
        std::cerr << "\n\ncqf" << c << "\n";
        cqf.dump_metadata();
        std::cerr << "\n\n";

        k = cqf.keybits() / 2;
        logger->info("Done loading cdbg. k is {}", k);
        logger->info("Iterating over cqf & building edgeSet ...");
        // build color class edges in a multi-threaded manner
        std::vector<std::thread> threads;
        for (uint32_t i = 0; i < nThreads; ++i) {
            threads.emplace_back(std::thread(&MSTMerger::buildPairedColorIdEdgesInParallel, this, i,
                                             std::ref(cqf), std::ref(tmpEdges), std::ref(curFileIdx),
                                             std::ref(cnts[i]),
                                             std::ref(maxId), std::ref(numOfKmers)));
        }
        for (auto &t : threads) { t.join(); }
    }
    if (not tmpEdges.empty()) {
        std::cerr << "tmpEdges.size(): " << tmpEdges.size() << " sizeof(element): " << sizeof(std::remove_reference<decltype(tmpEdges)>::type::value_type) <<  "\n";
        edgePairSortUniq(tmpEdges);
        std::cerr << "tmpEdges.size(): " << tmpEdges.size() << "\n";
        std::string filename(prefix+"tmp"+std::to_string(curFileIdx));
        std::cerr << "Filename: " << filename << "\n";
        std::ofstream tmpfile(filename, std::ios::out | std::ios::binary);
        uint64_t vecSize = tmpEdges.size();
        std::cerr << sizeof(vecSize) << "\n";
        tmpfile.write(reinterpret_cast<const char *>(&vecSize), sizeof(vecSize));
        std::cerr << "Done writing size\n";
        tmpfile.write(reinterpret_cast<const char *>(tmpEdges.data()),
                      sizeof(std::remove_reference<decltype(tmpEdges)>::type::value_type)*tmpEdges.size());
        tmpfile.close();
        curFileIdx++;
    }
    logger->info("Done writing the last tmp file");
    uint64_t totalEdges{0};
    for (uint32_t i = 0; i < nThreads; ++i) {
        totalEdges += cnts[i];
    }
    logger->info("Total number of kmers observed: {}, total number of edges (including duplicates): {}", numOfKmers, totalEdges);
    num_colorClasses = maxId + 1 + 1; // last one is for zer as the dummy color class with ID equal to actual num of color classes

    /*logger->info("Merge edges of different temp buckets in a sorted list and calculate weight on the fly.");
    // Add an edge between each color class ID and node zero
    logger->info("Adding edges from dummy node zero to each color class Id for {} color classes",
                 num_colorClasses);
    zero = static_cast<colorIdType>(num_colorClasses);
    for (colorIdType colorId = 0; colorId < num_colorClasses; colorId++) {
        edges->emplace_back(colorId, zero);
    }
    num_colorClasses++; // zero is now a dummy color class with ID equal to actual num of color classes
*/
    return true;
}

void MSTMerger::buildPairedColorIdEdgesInParallel(uint32_t threadId,
                                            CQF<KeyObject> &cqf,
                                            std::vector<std::pair<uint64_t, uint64_t>> &tmpEdges,
                                            uint64_t &curFileIdx,
                                            std::uint64_t& cnt,
                                            uint64_t &maxId, uint64_t &numOfKmers) {
    //std::cout << "THREAD ..... " << threadId << " " << cqf.range() << "\n";
    uint64_t kmerCntr{0}, localMaxId{0};
    __uint128_t startPoint = threadId * (cqf.range() / (__uint128_t) nThreads);
    __uint128_t endPoint =
            threadId + 1 == nThreads ? cqf.range() + 1 : (threadId + 1) * (cqf.range() / (__uint128_t) nThreads);
    /*std::cerr << threadId << ": s" << (uint64_t) (startPoint/(__uint128_t)0xFFFFFFFFFFFFFFFF) << " "
              << "sr" << (uint64_t) (startPoint%(__uint128_t)0xFFFFFFFFFFFFFFFF) << " "
            << "e" << (uint64_t) (endPoint/(__uint128_t)0xFFFFFFFFFFFFFFFF) << " "
            << "er" << (uint64_t) (endPoint%(__uint128_t)0xFFFFFFFFFFFFFFFF) << "\n";*/
    auto tmpEdgeListSize = MAX_ALLOWED_TMP_EDGES / nThreads;
    std::vector<std::pair<uint64_t , uint64_t >> edgeList;
    edgeList.reserve(tmpEdgeListSize);
    auto appendStore = [&]() {
        edgePairSortUniq(edgeList);
        colorMutex.lock();
        uint64_t prevSize=tmpEdges.size();
        tmpEdges.insert(tmpEdges.end(), edgeList.begin(), edgeList.end());
        if (tmpEdges.size() >= MAX_ALLOWED_TMP_EDGES_IN_FILE) {
            edgePairSortUniq(tmpEdges);
            std::string filename(prefix+"tmp"+std::to_string(curFileIdx));
            std::ofstream tmpfile(filename, std::ios::out | std::ios::binary);
            tmpfile.write(reinterpret_cast<const char *>(tmpEdges.size()), sizeof(size_t));
            tmpfile.write(reinterpret_cast<const char *>(tmpEdges.data()),
                    sizeof(std::remove_reference<decltype(tmpEdges)>::type::value_type)*tmpEdges.size());
            tmpfile.close();
            tmpEdges.clear();
            curFileIdx++;
        }
        colorMutex.unlock();
        cnt+=edgeList.size();
        edgeList.clear();
    };
    auto it = cqf.setIteratorLimits(startPoint, endPoint);
    while (!it.reachedHashLimit()) {
        KeyObject keyObject = *it;
        uint64_t curEqId = keyObject.count - 1;
        //nodes[curEqId] = 1; // set the seen color class id bit
        localMaxId = curEqId > localMaxId ? curEqId : localMaxId;
        // Add an edge between the color class and each of its neighbors' colors in dbg
        findNeighborEdges(cqf, keyObject, edgeList);
        if (edgeList.size() >= tmpEdgeListSize) {
            appendStore();
        }
        ++it;
        kmerCntr++;
        if (kmerCntr % 10000000 == 0) {
            std::cerr << "\rthread " << threadId << ": Observed " << (numOfKmers + kmerCntr) / 1000000 << "M kmers and " << cnt << " edges";
        }
    }
    appendStore();
    colorMutex.lock();
    maxId = localMaxId > maxId ? localMaxId : maxId;
    numOfKmers += kmerCntr;
    std::cerr << "\r";
    logger->info("Thread {}: Observed {} kmers and {} edges", threadId, numOfKmers, cnt);
    colorMutex.unlock();
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

    uint64_t nodeCnt1 = MSTQuery::getNodeCount(prefix1);
    uint64_t nodeCnt2 = MSTQuery::getNodeCount(prefix2);
    uint64_t mst1Zero{nodeCnt1 - 1}, mst2Zero{nodeCnt2 - 1};
    auto c1len{static_cast<uint64_t >(ceil(log2(nodeCnt1)))}, c2len{static_cast<uint64_t >(ceil(log2(nodeCnt2)))};
    uint64_t c2mask = (1ULL << c2len) - 1;
    std::cerr << "c1len: " << c1len << " c2len:" << c2len << " c2mask:" << c2mask << "\n";
    nonstd::optional<uint64_t> dummy{nonstd::nullopt};

    QueryStats dummyStats1, dummyStats2;

    logger->info("loaded the two msts with k={}. MST sizes are {}, {} respectively.", k, nodeCnt1,
                 nodeCnt2);
    std::ifstream cp(prefix + "newID2oldIDs");
    uint64_t cnt, cIdx;
    colorIdType c1, c2;
    cp.read(reinterpret_cast<char *>(&cnt), sizeof(cnt));
    logger->info("# of color classes based on count of colorPairs: {}", cnt);
    logger->info("# of color classes for mantis 1 : {} and for mantis 2 : {}", mst1Zero, mst2Zero);
    // colorPairs colors for mantis 1 and mantis 2 are 0-based --> color ID starts from 0
    colorPairs = sdsl::int_vector<>(cnt, 0, c1len+c2len);
    uint64_t maxIndex = 0;
    std::cerr << "Reading colorPairs from file " << cnt << "\n";
    for (auto i = 0; i < cnt; i++) {
        cp.read(reinterpret_cast<char *>(&cIdx), sizeof(cIdx));
        cp.read(reinterpret_cast<char *>(&c1), sizeof(c1));
        cp.read(reinterpret_cast<char *>(&c2), sizeof(c2));
        c1 = c1 == 0 ? mst1Zero : c1 - 1;
        c2 = c2 == 0 ? mst2Zero : c2 - 1;
        colorPairs[cIdx] = ((c1 << c2len) | c2);
        std::cerr << cIdx << " " << c1 << " " << c2 << " " << colorPairs[cIdx] << " " << (colorPairs[cIdx] >> c2len) << " " << (colorPairs[cIdx] & c2mask) << "\n";
        maxIndex = maxIndex>=cIdx?maxIndex:cIdx;
    }
    cp.close();


    mst1 = std::make_unique<MSTQuery>(prefix1, k, k, numOfFirstMantisSamples, logger);
    std::vector<colorIdType> colorsInCache;
    planCaching(mst1.get(), colorsInCache);
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

    mst2 = std::make_unique<MSTQuery>(prefix2, k, k, secondMantisSamples, logger);
    planCaching(mst2.get(), colorsInCache);
    logger->info("fixed cache size for mst2 is : {}", colorsInCache.size());
    //fillout fixed_cache2
    for (int64_t idx = colorsInCache.size() - 1; idx >= 0; idx--) {
        auto setbits = mst2->buildColor(colorsInCache[idx], dummyStats2, &lru_cache2[0], nullptr, &fixed_cache2, dummy);
        fixed_cache2[colorsInCache[idx]] = setbits;
    }
    logger->info("Done filling the fixed cache for mst2. Calling multi-threaded MSTBasedHammingDist .. ");

    uint64_t cntr{0};
    logger->info("MST 1 and 2 zeros: {}, {}", mst1Zero, mst2Zero);

    sdsl::int_vector<> ccSetBitCnts1(mst1Zero+1, 0, ceil(log2(numOfFirstMantisSamples)));
    sdsl::int_vector<> ccSetBitCnts2(mst2Zero+1, 0, ceil(log2(secondMantisSamples)));
    uint32_t maxWeightInFile = std::min(static_cast<uint32_t >(1000), numSamples);
    std::vector<uint64_t> bufferEndIndices;
    uint64_t size = splitHarmonically(MAX_ALLOWED_TMP_EDGES_IN_FILE/2, maxWeightInFile, bufferEndIndices);
    std::vector<uint64_t> bufferCurrIndices(bufferEndIndices.size(), 0);
    std::vector<uint64_t > bufferCnt(bufferEndIndices.size(), 0);
    for (auto i = 1; i < bufferCurrIndices.size(); i++) {
        bufferCurrIndices[i] = bufferEndIndices[i-1];
    }
    std::vector<std::pair<uint64_t , uint64_t >> output_buffer(size);
    std::vector<std::ofstream> outWeightFile;
    outWeightFile.reserve(maxWeightInFile);
    for (auto tmpOutIdx = 1; tmpOutIdx <= maxWeightInFile; tmpOutIdx++) {
        std::string filename = prefix+"w"+std::to_string(tmpOutIdx);
        outWeightFile.emplace_back(filename, std::ios::out | std::ios::binary);
        uint64_t cnt = 0;
        outWeightFile.back().write(reinterpret_cast<char*>(&cnt), sizeof(cnt));
    }
    std::vector<std::tuple<uint64_t , uint64_t , uint32_t >> inMemEdgeWeight;
    inMemEdgeWeight.reserve(1000000);
    auto store_edge_weight = [&outWeightFile,
                              &output_buffer,
                              &bufferEndIndices,
                              &bufferCurrIndices,
                              &bufferCnt,
                              &inMemEdgeWeight,
                              &maxWeightInFile](auto &pair, auto w) {
        w--;
        std::cerr << "w:" << w << " maxWeightInFile:" << maxWeightInFile << "\n";
        if (w < maxWeightInFile) {
            std::cerr << "curvsend: " << bufferCurrIndices[w] << " " << bufferEndIndices[w] << "\n";
            if (bufferCurrIndices[w] < bufferEndIndices[w]) {
                output_buffer[bufferCurrIndices[w]] = pair;
                bufferCurrIndices[w]++;
            } else {
                auto startIdx = w == 0?0:bufferEndIndices[w-1];
                outWeightFile[w].write(reinterpret_cast<char *>(output_buffer.data()+startIdx),
                                       sizeof(std::remove_reference<decltype(output_buffer)>::type::value_type) * (bufferEndIndices[w]-startIdx));
                bufferCnt[w] += (bufferEndIndices[w]-startIdx);
                bufferCurrIndices[w] = startIdx;
            }
        } else {
            inMemEdgeWeight.emplace_back(pair.first, pair.second, w);
        }
    };

    std::vector<TmpFileIterator> tis;
    Minheap_edge minheap;
    std::cerr << "curfileidx: " << curFileIdx << "\n";
    for (uint64_t i = 0; i < curFileIdx; i++) {
        tis.emplace_back(prefix+"tmp"+std::to_string(i), MAX_ALLOWED_TMP_EDGES_IN_FILE/(2.0*curFileIdx));
    }
    for (uint64_t i = 0; i < curFileIdx; i++) {
        if (tis[i].end()) continue;
        minheap.push(&tis[i]);
        std::cerr <<"iterator " << i << " pushed\n";
    }
    std::cerr << tis[0].get_val().first << " " << tis[0].get_val().second << "\n";
    std::cerr << minheap.empty() << "\n";
    std::cerr << minheap.top() << "\n";
    std::cerr << &tis[0] <<"\n";
    auto val = minheap.top()->get_val();
    while (!minheap.empty()) {
        do {
            auto cur = minheap.top();
            val = cur->get_val();
            std::cerr << "while " <<val.first << " " << val.second << " ";
            if (cur->next())
                minheap.replace_top(cur);
            else
                minheap.pop();
        } while (!minheap.empty() && val == minheap.top()->get_val());
        // calculate weight
        auto src = colorPairs[val.first];
        auto dest = colorPairs[val.second];
        auto src1 = src >> c2len;
        auto dest1 = dest >> c2len;
        auto w1 = mstBasedHammingDist(src1, dest1,
                                      mst1.get(), lru_cache1[0], queryStats1[0],
                                      fixed_cache1,
                                      ccSetBitCnts1);
        auto src2 = src & c2mask;
        auto dest2 = dest >> c2len;
        auto w2 = mstBasedHammingDist(src2, dest2,
                                      mst2.get(), lru_cache2[0], queryStats2[0],
                                      fixed_cache2,
                                      ccSetBitCnts2);
        auto w = w1 + w2;
        std::cerr << val.first << " " << val.second << " " << src << "->" << dest << " e1: " << src1 << "->" << dest1 << " e2: " << src2 << "->" << dest2 << " w: " << w1 << " " << w2 << " " << w << "\n";
        store_edge_weight(val, w);
    }
    for (auto cc : colorPairs) {
        auto cc1 = cc >> c2len;
        auto cc2 = cc & c2mask;
        auto w = ccSetBitCnts1[cc1] + ccSetBitCnts2[cc2];
        auto pair = std::make_pair(cc, zero);
        store_edge_weight(pair, w);
    }
    for (auto w = 0; w < outWeightFile.size(); w++) {
        auto startIdx = w == 0?0:bufferEndIndices[w-1];
        std::cerr << " bufferCurrIndices[" << w << "]:" << bufferCurrIndices[w] << " vs startIdx:" << startIdx << "\n";
        if (bufferCurrIndices[w] != startIdx) {
            outWeightFile[w].write(reinterpret_cast<char *>(output_buffer.data()+startIdx),
                                   sizeof(std::remove_reference<decltype(output_buffer)>::type::value_type) * (bufferCurrIndices[w]-startIdx));
            bufferCnt[w] += (bufferEndIndices[w]-startIdx);
        }
        std::cerr << "before writing to " << w << " cnt: " << bufferCnt[w] << "\n";
        outWeightFile[w].seekp(0);
        outWeightFile[w].write(reinterpret_cast<char*>(bufferCnt.data()+w), sizeof(uint64_t));
        outWeightFile[w].close();
    }

//    for (auto &edge : (*edges)) {
//        auto w1 = findWeight(edge, edge1list, srcEndIdx1, mst1Zero, true);
//        auto w2 = findWeight(edge, edge2list, srcEndIdx2, mst2Zero, false);
//        weightBuckets[w1 + w2 - 1]->push_back(edge);
//        if (++cntr % 10000000 == 0)
//            std::cerr << "\r" << cntr << " out of " << edges->size();
//    }
//    std::cerr << "\r";
//    edges.reset(nullptr);
    logger->info("Calculated the weight for the edges");
    return true;
}

/**
 * Finds Minimum Spanning Forest of color graph using Kruskal Algorithm
 *
 * The algorithm's basic implementation taken from
 * https://www.geeksforgeeks.org/kruskals-minimum-spanning-tree-using-stl-in-c/
 * @return List of connected components in the Minimum Spanning Forest
 */
void MSTMerger::kruskalMSF(AdjList * adjListPtr) {
    uint32_t bucketCnt = numSamples;
//    std::ofstream mstAdj(prefix + mantis::TEMP_MST_ADJ_FILE);
//    std::make_unique<std::vector<std::vector<std::pair<colorIdType, uint32_t> >>>(num_colorClasses);
    // Create disjoint sets
    DisjointTrees ds(num_colorClasses);
    uint64_t edgeCntr{0}, selectedEdgeCntr{0}, u, v, w;
//    uint32_t w{0};

    // Iterate through all sorted edges
    for (uint32_t bucketCntr = 0; bucketCntr < bucketCnt; bucketCntr++) {
        uint32_t edgeIdxInBucket = 0;
        w = bucketCntr + 1;
        TmpFileIterator it(prefix+"w"+std::to_string(w), MAX_ALLOWED_TMP_EDGES_IN_FILE);
//        for (auto &it : *(weightBuckets[bucketCntr])) {
        while (not it.end()) {
            auto edge = it.get_val();
            u = edge.first;
            v = edge.second;
            auto root_of_u = ds.find(u);
            auto root_of_v = ds.find(v);

            // Check if the selected edge is causing a cycle or not
            // (A cycle is induced if u and v belong to the same set)
            if (root_of_u != root_of_v) {
                // Merge two sets
                ds.merge(root_of_u, root_of_v, w);
                // Current edge will be in the MST
                adjListPtr->storeEdge(u, v, w);
//                mstAdj << u << " " << v << " " << w << "\n";
                mstTotalWeight += w;
                selectedEdgeCntr++;
            }
            edgeCntr++;
            if (edgeCntr % 1000000 == 0) {
                std::cerr << "\r" << edgeCntr << " edges processed and "
                          << selectedEdgeCntr << " were selected";
            }
            edgeIdxInBucket++;
            it.next();
            if (selectedEdgeCntr == num_colorClasses-1) break;
        }
        if (selectedEdgeCntr == num_colorClasses-1) break;
//        weightBuckets[bucketCntr].reset(nullptr);
    }
    std::cerr << "\r";
    mstTotalWeight++;//1 empty slot for root (zero)
    logger->info("MST Construction finished:"
                 "\n\t# of graph nodes: {}"
                 "\n\t# of graph edges: {}"
                 "\n\t# of merges (mst edges): {}"
                 "\n\tmst weight sum: {}",
                 num_colorClasses, edgeCntr, selectedEdgeCntr, mstTotalWeight);
//    mstAdj.close();
}

/**
 * calls kruskal algorithm to build an MST of the color graph
 * goes over the MST and fills in the int-vectors parentbv, bbv, and deltabv
 * serializes these three int-vectors as the encoding of color classes
 * @return true if encoding and serializing the DS is successful
 */
bool MSTMerger::encodeColorClassUsingMST() {
    // build mst of color class graph
    logger->info("before kruskal");
//    usleep(10000000);
    auto adjListPtr = std::make_unique<AdjList>(prefix, num_colorClasses, numSamples);
    kruskalMSF(adjListPtr.get());
    adjListPtr->loadCompactedAdjList();
    logger->info("after kruskal");
//    usleep(10000000);
    uint64_t nodeCntr{0};
    // encode the color classes using mst
    logger->info("Filling ParentBV...");
    sdsl::int_vector<> parentbv(num_colorClasses, 0, ceil(log2(num_colorClasses)));
    sdsl::bit_vector visited(num_colorClasses, 0);
    std::cerr << "after initializing parentbv and visited\n";
//    usleep(10000000);
    std::cerr << "zero: " << zero << " " << parentbv.size() << "\n";
    adjListPtr->hybridTreeWalk(zero, parentbv, visited);
    std::cerr << "\r";
    std::cerr << "Filled parentBV\n";
    uint64_t cntr = 0;
    for (auto i = 0; i < visited.size(); i++) {
        if (visited[i] == 0) {
            cntr++;
        }
    }
    std::cerr << "total " << visited.size() << ", visited: " << adjListPtr->visitedCnt << ", not visited: " << cntr << "\n";
    std::vector<uint64_t> thread_deltaOffset_and_parentEnd(nThreads, 0);
    uint64_t idx{0};
    uint64_t bucketSize = std::ceil(parentbv.size() / (double)nThreads);
    uint64_t doubleCheckTotWeight{0};
    for (auto i = 1; i <= adjListPtr->smallerSrcCnt.size(); i++) {
        auto limit = i == adjListPtr->smallerSrcCnt.size() ?
                adjListPtr->smallerSrcCnt.size() : adjListPtr->smallerSrcCnt[i];
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

    std::cerr << "\r";
    adjListPtr.reset(nullptr);
    std::cerr << "after deleting MST adjacency list\n";
    usleep(10000000);

    // fill in deltabv and bbv
    logger->info("Filling DeltaBV and BBV...");
    mst1.reset(new MSTQuery(prefix1, k, k, secondMantisSamples, logger));
    mst2.reset(new MSTQuery(prefix2, k, k, secondMantisSamples, logger));
//    mst1->loadIdx(prefix1);
//    mst2->loadIdx(prefix2);
    sdsl::bit_vector bbv(mstTotalWeight, 0);
    sdsl::int_vector<> deltabv(mstTotalWeight, 0, ceil(log2(numSamples)));
    std::cerr << "after initializing deltabv and bbv\n";
    usleep(10000000);

//    sdsl::bit_vector::select_1_type sbbv = sdsl::bit_vector::select_1_type(&bbv);
    std::vector<std::thread> threads;
    for (uint32_t t = 0; t < nThreads; ++t) {
        threads.emplace_back(std::thread(&MSTMerger::calcDeltasInParallel, this,
                                         t, thread_deltaOffset_and_parentEnd[t],
                                         std::ref(parentbv), std::ref(deltabv), std::ref(bbv), true));
    }
    for (auto &t : threads) { t.join(); }
    bbv[bbv.size()-1] = 1;
    std::cerr << "\r";

    std::cerr << "Done\n";
    usleep(10000000);

    logger->info("Serializing data structures parentbv, deltabv, & bbv...");
    sdsl::store_to_file(parentbv, std::string(prefix + mantis::PARENTBV_FILE));
    sdsl::store_to_file(deltabv, std::string(prefix + mantis::DELTABV_FILE));
    sdsl::store_to_file(bbv, std::string(prefix + mantis::BOUNDARYBV_FILE));
    logger->info("Done Serializing.");
    return true;
}

void MSTMerger::calcDeltasInParallel(uint32_t threadID, uint64_t deltaOffset,
                                     sdsl::int_vector<> &parentbv, sdsl::int_vector<> &deltabv,
                                     sdsl::bit_vector &bbv,
                                     bool isMSTBased) {

    uint64_t deltasKeptInMem = 1000;
    struct Delta {
        uint64_t startingOffset{0};
        std::vector<uint32_t> deltaVals;

        Delta() = default;

        Delta(uint64_t so) {
            startingOffset = so;
        }
    };
    std::vector<Delta> deltas;
    deltas.reserve(deltasKeptInMem);

    colorIdType mst1Zero = mst1->parentbv.size() - 1, mst2Zero = mst2->parentbv.size() - 1;
    auto c1len{static_cast<uint64_t >(ceil(log2(mst1Zero+1)))}, c2len{static_cast<uint64_t >(ceil(log2(mst2Zero+1)))};
    uint64_t c2mask = (1ULL << c2len) - 1;
    uint64_t bucketSize = std::ceil(parentbv.size() / (double)nThreads);
    uint64_t s = bucketSize * threadID;;
    uint64_t e = std::min(parentbv.size(), bucketSize * (threadID+1));
    for (uint64_t childIdx = s; childIdx < e; childIdx++) {
//        auto deltaOffset = (childIdx > 0) ? (sbbv(childIdx) + 1) : 0;
        deltas.emplace_back(deltaOffset);
        uint64_t child1{mst1Zero}, parent1{mst1Zero}, child2{mst2Zero}, parent2{mst2Zero};
        if (childIdx != zero) {
            auto child = colorPairs[childIdx];
            child1 = child >> c2len;
            child2 = child & c2mask;
        }
        if (parentbv[childIdx] != zero) {
            auto parent = colorPairs[parentbv[childIdx]];
            parent1 = parent >> c2len;
            parent2 = parent & c2mask;
        }

        auto firstDelta = getMSTBasedDeltaList(child1, parent1, mst1.get(), fixed_cache1,
                                               lru_cache1[threadID], queryStats1[threadID]);
        auto secondDelta = getMSTBasedDeltaList(child2, parent2, mst2.get(), fixed_cache2,
                                                lru_cache2[threadID], queryStats2[threadID]);
        deltas.back().deltaVals = firstDelta;
        for (auto &v : secondDelta) {
            v += numOfFirstMantisSamples;
            deltas.back().deltaVals.push_back(v);
        }
        deltaOffset += deltas.back().deltaVals.size();
        if (deltas.size() >= deltasKeptInMem) {
            colorMutex.lock();
            for (auto &v : deltas) {
                if (v.startingOffset > 0) {
                    bbv[v.startingOffset - 1] = 1;
                }
                for (auto cntr = 0; cntr < v.deltaVals.size(); cntr++) {
                    deltabv[v.startingOffset + cntr] = v.deltaVals[cntr];
                }
            }
            colorMutex.unlock();
            deltas.clear();
        }
    }

    colorMutex.lock();
    for (auto &v : deltas) {
        if (v.startingOffset > 0) {
            bbv[v.startingOffset - 1] = 1;
        }
        for (auto cntr = 0; cntr < v.deltaVals.size(); cntr++) {
            deltabv[v.startingOffset + cntr] = v.deltaVals[cntr];
        }
    }
    colorMutex.unlock();
}

/**
 * finds the neighbors of each kmer in the cqf,
 * and adds an edge of the element's colorId and its neighbor's (0-based)
 * @param cqf (required to query for existence of neighbors)
 * @param it iterator to the elements of cqf
 */
void MSTMerger::findNeighborEdges(CQF<KeyObject> &cqf, KeyObject &keyobj, std::vector<std::pair<uint64_t , uint64_t> > &edgeList) {
    dna::canonical_kmer curr_node(static_cast<int>(k), keyobj.key);
    workItem cur = {curr_node, static_cast<colorIdType>(keyobj.count - 1)};
    uint64_t neighborCnt{0};
    for (auto &nei : neighbors(cqf, cur)) {
        neighborCnt++;
        if (cur.colorId < nei.colorId) {
//            Edge e(static_cast<colorIdType>(cur.colorId), static_cast<colorIdType>(nei.colorId));

            edgeList.emplace_back(cur.colorId, nei.colorId);
        }
    }
}

/**
 * finds neighbors of a node in cqf
 * @param cqf
 * @param n : work_item containing node and colorId (colorId will be filled)
 * @return set of neighbors for current node n and their colorIds (0-based)
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
    if ((not fixed_cache.empty()) and fixed_cache.find(eqid) != fixed_cache.end()) {
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
                                        QueryStats &queryStats,
                                        std::unordered_map<uint64_t, std::vector<uint64_t>> &fixed_cache,
                                        sdsl::int_vector<> &ccSetBitCnt) {

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

std::vector<uint32_t> MSTMerger::getMSTBasedDeltaList(uint64_t eqid1, uint64_t eqid2,
                                                      MSTQuery * mstPtr,
                                                      std::unordered_map<uint64_t, std::vector<uint64_t>>& fixed_cache,
                                                      LRUCacheMap &lru_cache,
                                                      QueryStats &queryStats) {
    std::vector<uint32_t> res;
    if (eqid1 == eqid2) return res;
    std::vector<uint64_t> eq1, eq2;
    buildMSTBasedColor(eqid1, mstPtr, lru_cache, eq1, queryStats, fixed_cache);
    buildMSTBasedColor(eqid2, mstPtr, lru_cache, eq2, queryStats, fixed_cache);
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

void MSTMerger::planCaching(MSTQuery *mst,std::vector<colorIdType> &colorsInCache) {
    std::vector<Cost> mstCost(mst->parentbv.size());

    logger->info("In planner ..");
    // setting local edge costs
    // local cost for each node is its degree (in + out)
    uint64_t src{0}, edgeCntr{0};
    uint64_t queryCost = mst->parentbv.size()-1;
    for (auto i = 0; i < queryCost; i++) {
        mstCost[i].numQueries = queryCost;
        queryCost = queryCost/(i+1);
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
