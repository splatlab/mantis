/*
 * ============================================================================
 *
 *         Authors:
 *         Fatemeh Amodaresi, falmodar@cs.umd.edu
 *         Jamshed Khan, jamshed@cs.umd.edu
 *   Organization:  University of Maryland
 *
 * ============================================================================
 */

#include "cqfMerger.h"
#include "tbb/parallel_sort.h"

template<typename qf_obj, typename key_obj>
CQF_merger<qf_obj, key_obj>::
CQF_merger(std::string &firstCQFdir, std::string &secondCQFdir,
            std::string &outputCQFdir, spdlog::logger *logger, uint64_t threadNum)
            : console(logger), threadCount(threadNum)
{
    console -> info("Loading metadata for the first input colored dBG from disk.");

    cdbg1 = ColoredDbg<SampleObject<CQF<KeyObject> *>, KeyObject>(firstCQFdir, MANTIS_DBG_IN_MEMORY);//MANTIS_DBG_ON_DISK);

    console -> info("Read colored dBG over {} samples, with {} cqf files and {} color-class files.",
                    cdbg1.get_num_samples(), cdbg1.get_numBlocks(), cdbg1.get_eq_class_file_count());


    console -> info("Loading metadata for the second input colored dBG from disk.");

    cdbg2 = ColoredDbg<SampleObject<CQF<KeyObject> *>, KeyObject>(secondCQFdir, MANTIS_DBG_IN_MEMORY);//MANTIS_DBG_ON_DISK);

    console -> info("Read colored dBG over {} samples, with {} cqf files and {} color-class files.",
                    cdbg2.get_num_samples(), cdbg2.get_numBlocks(), cdbg2.get_eq_class_file_count());


    if (cdbg1.get_current_cqf() == nullptr) {
        console -> error("The first blocked cqf for first input mantis is null.");
        std::exit(3);
    }
    if (cdbg2.get_current_cqf() == nullptr) {
        console -> error("The first blocked cqf for second input mantis is null.");
        std::exit(3);
    }
    if(!cdbg1.get_current_cqf() -> check_similarity(cdbg2.get_current_cqf()))
    {
        console -> error("The CQF files of the colored dBGs are not similar.");
        exit(1);
    }

    console -> info("Initializing the output Mantis.");
    cdbg = ColoredDbg<SampleObject<CQF<KeyObject> *>, KeyObject>(cdbg1, cdbg2, outputCQFdir, MANTIS_DBG_IN_MEMORY);//MANTIS_DBG_ON_DISK);
    colorBits[0] = ceil(log2(cdbg1.get_num_eqclasses()));
    colorBits[1] = ceil(log2(cdbg2.get_num_eqclasses()));
    console->info("colorCnt1: {}, colorBits1:{} || colorCnt2: {}, colorBits2: {}", cdbg1.get_num_eqclasses(), colorBits[0], cdbg2.get_num_eqclasses(), colorBits[1]);
    start_time_ = std::time(nullptr);
    kbits = cdbg1.get_current_cqf()->keybits();
    hashmode = cdbg1.get_current_cqf()->hash_mode();
    seed = cdbg1.get_current_cqf()->seed();
    kmerMask = (1ULL << kbits) - 1;

//    minimizerKeyList[0].resize(cdbg1.minimizerBlock.size());
//    minimizerKeyList[1].resize(cdbg2.minimizerBlock.size());
//    minimizerColorList[0].resize(cdbg1.minimizerBlock.size());
//    minimizerColorList[1].resize(cdbg2.minimizerBlock.size());
//    for (auto &m : minimizerKeyColorList[0]) m.reset(new std::vector<std::pair<uint64_t, colorIdType>>());
//    for (auto &m : minimizerKeyColorList[1]) m.reset(new std::vector<std::pair<uint64_t, colorIdType>>());
}

template <typename qf_obj, typename key_obj>
uint64_t CQF_merger<qf_obj, key_obj>::
walkBlockedCQF(ColoredDbg<qf_obj, key_obj> &curCdbg, const uint64_t curBlock, bool isSecond) {
    uint64_t minMinimizer{invalid}, maxMinimizer{MAX_MINIMIZER};
    curCdbg.replaceCQFInMemory(curBlock);

    if (curBlock < curCdbg.get_numBlocks()) {
        const CQF<key_obj> *cqf1 = curCdbg.get_current_cqf();
        cqf1->dump_metadata();
        uint64_t count{cqf1->dist_elts()};
        auto minmax = curCdbg.getMinMaxMinimizer(curBlock);
        minMinimizer = minmax.first;
        maxMinimizer = minmax.second;
        std::vector<std::vector<uint32_t>> minimizerIdx(threadCount);
        for (auto &v: minimizerIdx)
            v.resize(maxMinimizer-minMinimizer+1, 0);
        std::vector<std::thread> threads;
        if (threadCount > 1) {
            for (auto threadId = 0; threadId < threadCount; threadId++) {
                threads.emplace_back(std::thread([&, threadId]() {
                    __uint128_t startPoint = threadId * (cqf1->range() / (__uint128_t) threadCount);
                    __uint128_t endPoint =
                            threadId + 1 == threadCount ? cqf1->range() + 1 : (threadId + 1) * (cqf1->range() /
                                                                                                (__uint128_t) threadCount);
                    auto it = cqf1->setIteratorLimits(startPoint, endPoint);
                    while (!it.reachedHashLimit()) {
                        auto keyval = it.get_cur_hash();
                        auto key = hash_64i(keyval.key, kmerMask);
                        auto minimizerPair = curCdbg.findMinimizer(key, kbits);
                        if (minimizerPair.first == minimizerPair.second) {
                            minimizerPair.second = invalid;
                        }
                        std::vector<uint64_t> pairs{minimizerPair.first, minimizerPair.second};
                        for (auto minimizer : pairs) {
                            if (minimizer != invalid and minimizer >= minMinimizer and minimizer <= maxMinimizer) {
                                minimizerIdx[threadId][minimizer - minMinimizer]++;
                            }
                        }
                        ++it;
                    }
                }));
            }
            for (auto &t : threads) { t.join(); }
            threads.clear();
            for (auto m = 0; m <= maxMinimizer - minMinimizer; m++) {
                uint32_t prevCnt = 0, cur = 0;
                for (auto idx = 0; idx < threadCount; idx++) {
                    cur = minimizerIdx[idx][m];
                    minimizerIdx[idx][m] = prevCnt;
                    prevCnt += cur;
                }
            }
        }
        uint64_t end = curCdbg.minimizerBlock.size()-1;
        for (auto i = minMinimizer; i <= maxMinimizer; i++) {
            minimizerKeyList[isSecond][end-i] = std::make_unique<std::vector<uint64_t>>(curCdbg.minimizerCntr[i], 0);
            minimizerColorList[isSecond][end-i] = std::make_unique<sdsl::int_vector<>>(curCdbg.minimizerCntr[i], 0, colorBits[isSecond]);
        }
//        console->info("Start walking Cqf{}. isSecond:{}", curBlock, isSecond);
        std::vector<SpinLockT> mutexes(maxMinimizer-minMinimizer+1);

        for (auto threadId = 0; threadId < threadCount; threadId++) {
            threads.emplace_back(std::thread([&, threadId](){
                __uint128_t startPoint = threadId * (cqf1->range() / (__uint128_t) threadCount);
                __uint128_t endPoint =
                        threadId + 1 == threadCount ? cqf1->range() + 1 : (threadId + 1) * (cqf1->range() / (__uint128_t) threadCount);
//                std::stringstream ss;
//                ss << "threadId: " << std::to_string(threadId) << "\n";
//                std::cerr << ss.str();
                auto it = cqf1->setIteratorLimits(startPoint, endPoint);
                uint64_t tmpSize{100000}, cntr{0};
                std::vector<KeyColorMin> tmpList;
                tmpList.reserve(tmpSize);
                while (!it.reachedHashLimit()) {
                    auto keyval = it.get_cur_hash();
                    auto key = hash_64i(keyval.key, kmerMask);
                    auto minimizerPair = curCdbg.findMinimizer(key, kbits);
                    if (minimizerPair.first == minimizerPair.second) {
                        minimizerPair.second = invalid;
                    }
                    std::vector<uint64_t> pairs{minimizerPair.first, minimizerPair.second};
                    for (auto minimizer : pairs) {
                        if (minimizer != invalid and minimizer >= minMinimizer and minimizer <= maxMinimizer) {
                            tmpList.emplace_back(keyval.key, keyval.count, minimizer);
                            cntr++;
                        }
                    }
                    ++it;
                    if (tmpList.size() >= tmpSize or it.reachedHashLimit()) {
                        tbb::parallel_sort(tmpList.begin(), tmpList.end(),
                                      [](const KeyColorMin &f, const KeyColorMin &s) {
                                          return f.minimizer == s.minimizer? f.key < s.key : f.minimizer < s.minimizer;
                                      });
//                        std::sort(std::execution::par_unseq, tmpList.begin(), tmpList.end(),
//                                [](const KeyColorMin &f, const KeyColorMin &s) {
//                            return f.minimizer == s.minimizer? f.key < s.key : f.minimizer < s.minimizer;
//                        });
                        uint64_t i = 0;
                        while (i < tmpList.size()) {
                            auto prev = tmpList[i].minimizer;
                            mutexes[prev - minMinimizer].lock();
                            while (i < tmpList.size() and tmpList[i].minimizer == prev) {
                                int64_t curIdx = minimizerIdx[threadId][prev - minMinimizer];
                                (*minimizerKeyList[isSecond][end-prev])[curIdx] = tmpList[i].key;
                                (*minimizerColorList[isSecond][end-prev])[curIdx] = tmpList[i].color;
                                minimizerIdx[threadId][prev - minMinimizer]++;
                                i++;
                            }
                            mutexes[prev - minMinimizer].unlock();
                        }
                        tmpList.clear();
                    }
                }
            }));
        }
        for (auto &t : threads) { t.join(); }
//        console->info("Done walking Cqf{}. isSecond:{}", curBlock, isSecond);
//        std::cerr << "\r" << (isSecond?"Second ":"First ") << cntr << " (main & duplicate) inserts for " << count << " kmers\n";
        // validating the minimizer observed count;
        /*for (auto i = minMinimizer; i <= maxMinimizer; i++) {
            if (minimizerIdx[numThreads-1][i - minMinimizer] != curCdbg.minimizerCntr[i]) {
                console->error("Did not observe all the kmers with minimizer {}. Observed:{}, Expected{}", i, minimizerIdx[numThreads-1][i-minMinimizer], curCdbg.minimizerCntr[i]);
                std::exit(3);
            }
            for (auto idx = 1; idx < minimizerKeyList[isSecond][i]->size(); idx++) {
                if ((*minimizerKeyList[isSecond][i])[idx-1] >= (*minimizerKeyList[isSecond][i])[idx]) {
                    std::cerr << "didn't work yet\n";
                    std::exit(3);
                }
            }
        }*/
    }
    curCdbg.replaceCQFInMemory(invalid);
    return maxMinimizer;
}

template <typename qf_obj, typename key_obj>
uint64_t CQF_merger<qf_obj, key_obj>::
sample_colorID_pairs(uint64_t sampleKmerCount)
{
    std::cerr << " sampleColorIDPairs qbits: " << cdbg.get_qbits() << "\n";

    auto t_start = time(nullptr);

    uint64_t sizes[2];
    sizes[0] = cdbg1.minimizerBlock.size();
    sizes[1] = cdbg2.minimizerBlock.size();
    for (auto i = 0; i < 2; i++) {
        minimizerKeyList[i] = std::vector<std::unique_ptr<std::vector<uint64_t >>>(sizes[i]);
        minimizerColorList[i] = std::vector<std::unique_ptr<sdsl::int_vector<>>>(sizes[i]);
    }

    console -> info("Sampling at least {} kmers. Time-stamp = {}.",
                    sampleKmerCount, time(nullptr) - start_time_);

    uint64_t kmerCount = 0;

    uint64_t curBlock{0};
    console->info("cdbg1.numBlocks={}, cdbg2.numBlocks={}", cdbg1.get_numBlocks(), cdbg2.get_numBlocks());
    uint64_t maxMinimizer{0}, minMinimizer{0},
            maxMinimizer1{0}, maxMinimizer2{0};
    while(kmerCount < sampleKmerCount and
          (curBlock < cdbg1.get_numBlocks() or curBlock < cdbg2.get_numBlocks())) {
        console->info("Current block={}", curBlock);
        maxMinimizer1 = walkBlockedCQF(cdbg1, curBlock, false);
        maxMinimizer2 = walkBlockedCQF(cdbg2, curBlock, true);

        maxMinimizer =
                // If it's the last block for both of the CQFs
                curBlock >= cdbg1.get_numBlocks() - 1 and curBlock >= cdbg2.get_numBlocks() - 1 ? std::max(maxMinimizer1, maxMinimizer2) :
                curBlock >= cdbg1.get_numBlocks() ? maxMinimizer2 : // If we've passed the last block for CQF1
                curBlock >= cdbg2.get_numBlocks() ? maxMinimizer1 : // If we've passed the last block for CQF2
                std::min(maxMinimizer1, maxMinimizer2);
        std::cerr << "\rMin minimizer=" << minMinimizer << " Max minimizer=" << maxMinimizer << "\n";
        curBlock++;
        uint64_t end[2];
        end[0] = sizes[0]-1;
        end[1] = sizes[1]-1;
        for (uint64_t b = minMinimizer; b <= maxMinimizer; b++) {
            // merge the two keys from cqf1 and cqf2
            uint64_t cntr0{0}, cntr1{0}, size0{0}, size1{0};
            if (minimizerKeyList[0][end[0]-b])
                size0 = minimizerKeyList[0][end[0]-b]->size();
            if (minimizerKeyList[1][end[1]-b])
                size1 = minimizerKeyList[1][end[1]-b]->size();
            if (b % 500 == 0)
                std::cerr << "\rminimizer " << b
                          << " size=" << size0 << "," << size1 << "     ";
            while (cntr0 < size0 and cntr1 < size1) {
                uint64_t key0 = (*minimizerKeyList[0][end[0]-b])[cntr0];
                uint64_t color0 = (*minimizerColorList[0][end[0]-b])[cntr0];
                uint64_t key1 = (*minimizerKeyList[1][end[1]-b])[cntr1];
                uint64_t color1 = (*minimizerColorList[1][end[1]-b])[cntr1];
                /*if (kmerCount < 20000) {
                    std::cerr << ">>false " << key0 << " " << cntr0 << "\n";
                    std::cerr << ">>true " << key1 << " " << cntr1 << "\n";
                }*/
                if (key0 < key1) {
                    sampledPairs[std::make_pair(color0, 0)]++;
                    ++cntr0;
                } else if (key0 > key1) {
                    sampledPairs[std::make_pair(0, color1)]++;
                    ++cntr1;
                } else { // it0->first == it1->first
                    sampledPairs[std::make_pair(color0, color1)]++;
                    ++cntr0;
                    ++cntr1;
                }
                cdbg.minimizerCntr[b]++;
                kmerCount++;
            }
            if (cntr0 == size0) {
                while (cntr1 < size1) {
                    uint64_t color = (*minimizerColorList[1][end[1]-b])[cntr1];
                    sampledPairs[std::make_pair(0, color)]++;
                    cdbg.minimizerCntr[b]++;
                    kmerCount++;
                    ++cntr1;
                }
            } else {
                while (cntr0 < size0) {
                    uint64_t color = (*minimizerColorList[0][end[0]-b])[cntr0];
                    sampledPairs[std::make_pair(color, 0)]++;
                    cdbg.minimizerCntr[b]++;
                    kmerCount++;
                    ++cntr0;
                }
            }
        }
        for (auto i = 0; i < 2; ++i) {
            std::cerr << minimizerKeyList[i].size();
            minimizerKeyList[i].erase(minimizerKeyList[i].begin()+(end[i]-maxMinimizer), minimizerKeyList[i].end());
            minimizerColorList[i].erase(minimizerColorList[i].begin()+(end[i]-maxMinimizer), minimizerColorList[i].end());
            std::cerr << " to " << minimizerKeyList[i].size() << " ";
        }
        minMinimizer = maxMinimizer+1;
        std::cerr << "\r";
    }
    /*for (auto i = 0; i < 2; i++) {
        minimizerKeyList[i].clear();
        minimizerColorList[i].clear();
        minimizerKeyList[i].shrink_to_fit();
        minimizerColorList[i].shrink_to_fit();
    }*/
    console -> info("Sampled {} k-mers, color-classes found: {}. Time-stamp = {}.",
                    kmerCount, sampledPairs.size(), time(nullptr) - start_time_);
    typedef std::pair<uint64_t, std::pair<uint64_t, uint64_t>> CountAndIdPair;
    // By default a max heap is created ordered
    // by first element of pair.
    std::priority_queue<CountAndIdPair, std::vector<CountAndIdPair>> maxPQ;
    for(auto p = sampledPairs.begin(); p != sampledPairs.end(); ++p) {
        maxPQ.push(std::make_pair(p->second, p->first));
    }

    uint64_t colorId = 0;
    while(!maxPQ.empty())
    {
        sampledPairs[maxPQ.top().second] = colorId;
        maxPQ.pop();
        colorId++;
    }

    console -> info("Sampled {} color-id pairs. Time-stamp = {}.", sampledPairs.size(),
                    time(nullptr) - start_time_);


    auto t_end = time(nullptr);
    console -> info("Sampling abundant color-id pairs took time {} seconds.", t_end - t_start);
    return curBlock;
}

template <typename qf_obj, typename key_obj>
uint64_t CQF_merger<qf_obj, key_obj>::
store_colorID_pairs(uint64_t startingBlock)
{
    if (startingBlock == 0) {
        console->error("Should have at least walked over one cqf block to sample kmers. But {} blocks have been processed.", startingBlock);
        std::exit(3);
    }
    auto t_start = time(nullptr);

    console -> info("Writing the non-sampled color-id pairs to disk-files of form ({}). Time-stamp = {}.",
                    EQ_ID_PAIRS_FILE + std::string("_X_Y"), time(nullptr) - start_time_);

    uint64_t writtenPairsCount = 0;

    std::string output = cdbg.prefix + EQ_ID_PAIRS_FILE;
    console -> info("Iterating over the CQFs for the non-sampled color-id pairs starting from block {}."
                    " writing in file \"{}\"", startingBlock, output);
    std::ofstream diskBucket(output, std::ios::out);
    uint64_t curBlock1{startingBlock}, curBlock2{startingBlock}, kmerCount{0};
    uint64_t maxMinimizer{0}, minMinimizer;
    minMinimizer = curBlock1 >= cdbg1.get_numBlocks() and curBlock2 >= cdbg2.get_numBlocks() ? 0 :
            curBlock1 >= cdbg1.get_numBlocks() ? cdbg2.getMinMaxMinimizer(curBlock2).first :
            curBlock2 >= cdbg2.get_numBlocks() ? cdbg1.getMinMaxMinimizer(curBlock1).first :
            std::min(cdbg1.getMinMaxMinimizer(curBlock1).first, cdbg2.getMinMaxMinimizer(curBlock2).first);

    uint64_t maxMinimizer1{cdbg1.getMinMaxMinimizer(curBlock1-1).second}, maxMinimizer2{cdbg2.getMinMaxMinimizer(curBlock2-1).second};
    cdbg.initializeNewCQFBlock(invalid, kbits, cdbg.get_qbits(), hashmode, seed);
    uint64_t ends[2];
    ends[0] = cdbg1.minimizerBlock.size()-1;
    ends[1] = cdbg2.minimizerBlock.size()-1;
    while(curBlock1 < cdbg1.get_numBlocks() or curBlock2 < cdbg2.get_numBlocks()) {
        console->info("Current block={},{}", curBlock1, curBlock2);
        if (curBlock2 >= cdbg2.get_numBlocks() or
                (curBlock1 < cdbg1.get_numBlocks() and maxMinimizer1 <= maxMinimizer2)) {
            maxMinimizer1 = walkBlockedCQF(cdbg1, curBlock1, false);
            curBlock1++;
        } else {
            maxMinimizer2 = walkBlockedCQF(cdbg2, curBlock2, true);
            curBlock2++;
        }
        maxMinimizer =
                // If it's the last block for both of the CQFs
                curBlock1 >= cdbg1.get_numBlocks() and curBlock2 >= cdbg2.get_numBlocks() ? std::max(maxMinimizer1, maxMinimizer2) :
                curBlock1 >= cdbg1.get_numBlocks() ? maxMinimizer2 : // If we've passed the last block for CQF1
                curBlock2 >= cdbg2.get_numBlocks() ? maxMinimizer1 : // If we've passed the last block for CQF2
                std::min(maxMinimizer1, maxMinimizer2);
        std::cerr << "\rMin minimizer=" << minMinimizer << " Max minimizer=" << maxMinimizer << "\n";
        uint64_t clearedBytes = 0;

        for (uint64_t b = minMinimizer; b <= maxMinimizer; b++) {
            // merge the two keys from cqf1 and cqf2
            uint64_t cntr0{0}, cntr1{0}, size0{0}, size1{0};
            if (minimizerKeyList[0][ends[0]-b])
                size0 = minimizerKeyList[0][ends[0]-b]->size();
            if (minimizerKeyList[1][ends[1]-b])
                size1 = minimizerKeyList[1][ends[1]-b]->size();
            if (b % 500 == 0)
                std::cerr << "\rminimizer " << b << " " << ends[0]-b << " " << ends[1]-b
                          << " size=" << size0 << "," << size1 << "     ";
            while (cntr0 < size0 and cntr1 < size1) {
                uint64_t key0 = (*minimizerKeyList[0][ends[0]-b])[cntr0];
                uint64_t color0 = (*minimizerColorList[0][ends[0]-b])[cntr0];
                uint64_t key1 = (*minimizerKeyList[1][ends[1]-b])[cntr1];
                uint64_t color1 = (*minimizerColorList[1][ends[1]-b])[cntr1];
                auto colorPair = std::make_pair(color0, color1);
                if (key0 < key1) {
                    colorPair = std::make_pair(color0, 0);
                    ++cntr0;
                } else if (key0 > key1) {
                    colorPair = std::make_pair(0, color1);
                    ++cntr1;
                } else { // keys are the same for the two CQFs : it0->first == it1->first
                    ++cntr0;
                    ++cntr1;
                }

                if(sampledPairs.find(colorPair) == sampledPairs.end())
                {
                    diskBucket << colorPair.first << " " << colorPair.second << "\n";
                    writtenPairsCount++;
                }
                cdbg.minimizerCntr[b]++;
                kmerCount++;
            }
            if (cntr0 == size0) {
                while (cntr1 < size1) {
                    uint64_t color = (*minimizerColorList[1][ends[1]-b])[cntr1];
                    auto colorPair = std::make_pair(0ULL, color);
                    if(sampledPairs.find(colorPair) == sampledPairs.end())
                    {
                        diskBucket << colorPair.first << " " << colorPair.second << "\n";
                        writtenPairsCount++;
                    }
                    cdbg.minimizerCntr[b]++;
                    kmerCount++;
                    ++cntr1;
                }
            } else {
                while (cntr0 < size0) {
                    uint64_t color = (*minimizerColorList[0][ends[0]-b])[cntr0];
                    auto colorPair = std::make_pair(color, 0ULL);
                    if(sampledPairs.find(colorPair) == sampledPairs.end())
                    {
                        diskBucket << colorPair.first << " " << colorPair.second << "\n";
                        writtenPairsCount++;
                    }
                    cdbg.minimizerCntr[b]++;
                    kmerCount++;
                    ++cntr0;
                }
            }
            clearedBytes += ((size0* (64+colorBits[0]))+(size1* (64+colorBits[1])));
            for (auto i = 0; i < 2; ++i) {
                minimizerKeyList[i].erase(minimizerKeyList[i].begin()+(ends[i]-b), minimizerKeyList[i].end());
                minimizerColorList[i].erase(minimizerColorList[i].begin()+(ends[i]-b), minimizerColorList[i].end());
                minimizerKeyList[i].shrink_to_fit();
                minimizerColorList[i].shrink_to_fit();
            }
        }
        /*for (auto i = 0; i < 2; i++) {
            minimizerKeyList[i].shrink_to_fit();
            minimizerColorList[i].shrink_to_fit();
        }*/
        std::cerr << "\r";
        std::cerr << "Memory cleared in this round: " << std::round((clearedBytes) / std::pow(2, 23)) << "MB\n";
        minMinimizer = maxMinimizer+1;
    }
    for (auto i = 0; i < 2; i++) {
        minimizerKeyList[i].clear();
        minimizerColorList[i].clear();
        minimizerKeyList[i].shrink_to_fit();
        minimizerColorList[i].shrink_to_fit();
    }
    console -> info("Distinct kmers found {}, color-id pairs written to disk {}. Time-stamp = {}.",
                    kmerCount, writtenPairsCount, time(nullptr) - start_time_);

    diskBucket.flush();
    diskBucket.close();

    // force currentCQFBlock to get loaded into memory
    std::cerr << "1.5. Before sortUnique colorPairs\n\n";
    auto t_end = time(nullptr);
    console -> info("Filling up the disk-buckets with color-id pairs took time {} seconds.", t_end - t_start);
    sortUniq_colorID_pairs();

    return kmerCount;
}

template <typename qf_obj, typename key_obj>
void CQF_merger<qf_obj, key_obj>::
sortUniq_colorID_pairs()
{
    if(!system(nullptr))
    {
        console -> error("Command processor does not exist.");
        exit(1);
    }
    uint64_t maxMemoryForSort = 10;

    std::string diskBucket = cdbg.prefix + EQ_ID_PAIRS_FILE;
    console -> info("Filtering out the unique eq-id pairs from file {} with {} threads. Time-stamp = {}",
                    diskBucket, threadCount, time(nullptr) - start_time_);
    std::string sysCommand = "sort -t' ' -u -n -k1,1 -k2,2";
    sysCommand += " --parallel=" + std::to_string(threadCount);
    sysCommand += " -S " + std::to_string(maxMemoryForSort) + "G";
    sysCommand += " -o " + diskBucket + " " + diskBucket;

    console -> info("System command used:\n{}", sysCommand);
    system(sysCommand.c_str());
}

template <typename qf_obj, typename key_obj>
void CQF_merger<qf_obj, key_obj>::
build_MPHF()
{
    auto t_start = time(nullptr);

    console -> info("Building an MPH (Minimal Perfect Hash) table per disk-bucket. Time-stamp = {}.",
                    time(nullptr) - start_time_);
    // Lowest bit/elem is achieved with gamma=1, higher values lead to larger mphf but faster construction/query.
    double gammaFactor = 2.0;	// gamma = 2 is a good tradeoff (leads to approx 3.7 bits/key).

    // Build the mphf.
    std::string colorIdPairFile = cdbg.prefix + EQ_ID_PAIRS_FILE;
    std::ifstream input(colorIdPairFile);
    // new lines will be skipped unless we stop it from happening:
    input.unsetf(std::ios_base::skipws);
    // count the newlines with an algorithm specialized for counting:
    auto colorIdPairCount = static_cast<unsigned int>(std::count(
                std::istream_iterator<char>(input),
                std::istream_iterator<char>(),
                '\n'));
    input.close();

    console -> info("Total non-sampled colorId pair count: {}", colorIdPairCount);
    if (colorIdPairCount > 0) {
        ColorIdPairIterator kb(colorIdPairFile);
        ColorIdPairIterator ke(colorIdPairFile, true);
        auto colorPairIt = boomphf::range(kb, ke);
        colorMph = std::make_unique<boophf_t>(colorIdPairCount, colorPairIt, threadCount, gammaFactor);
        console->info("Total memory consumed by all the MPH tables = {} MB.",
                      (colorMph->totalBitSize() / 8) / (1024 * 1024));
        auto t_end = time(nullptr);
        console->info("Building the MPH tables took time {} seconds and {} memory", t_end - t_start,
                      colorMph->totalBitSize());
    } else {
        console->info("No non-sampled colorIds");
    }
}

template <typename qf_obj, typename key_obj>
void CQF_merger<qf_obj, key_obj>::build_CQF()
{
    auto t_start = time(nullptr);
    uint64_t ends[2];
    ends[0] = cdbg1.minimizerBlock.size();
    ends[1] = cdbg2.minimizerBlock.size();
    for (auto i = 0; i < 2; i++) {
        minimizerKeyList[i] = std::vector<std::unique_ptr<std::vector<uint64_t >>>(ends[i]);
        minimizerColorList[i] = std::vector<std::unique_ptr<sdsl::int_vector<>>>(ends[i]);
        ends[i]--;
    }
    console -> info("At CQFs merging phase. Time-stamp = {}.\n", time(nullptr) - start_time_);

    std::vector<std::pair<uint64_t, uint64_t>> tmp_kmers;
    // reserve 1/3rd more than the threshold as the actual count of kmers is gonna be +-epsilon of the threshold
    tmp_kmers.reserve(cdbg.getCqfSlotCnt());

    uint64_t kmerCount{0}, foundAbundantId{0};
    uint64_t curBlock1{0}, curBlock2{0}, outputCQFBlockId{0}, occupiedSlotsCnt{0}, currMinimizerSlotsCnt{0};

    uint64_t maxMinimizer{0}, minMinimizer{0};
    uint64_t colorId{0}, epsilon{100}, maxMinimizerCnt{0};
    KeyObject keyObj;

    uint64_t maxMinimizer1{0}, maxMinimizer2{0};
    cdbg.initializeNewCQFBlock(invalid, kbits, cdbg.get_qbits(), hashmode, seed);
    while(curBlock1 < cdbg1.get_numBlocks() or curBlock2 < cdbg2.get_numBlocks()) {
        console->info("Current block={},{}", curBlock1, curBlock2);
        if (curBlock1 == 0 and curBlock2 == 0) {
            maxMinimizer1 = walkBlockedCQF(cdbg1, curBlock1, false);
            maxMinimizer2 = walkBlockedCQF(cdbg2, curBlock2, true);
            curBlock1++;
            curBlock2++;
        }
        else if (curBlock2 >= cdbg2.get_numBlocks() or
            (curBlock1 < cdbg1.get_numBlocks() and maxMinimizer1 <= maxMinimizer2)) {
            maxMinimizer1 = walkBlockedCQF(cdbg1, curBlock1, false);
            curBlock1++;
        } else {
            maxMinimizer2 = walkBlockedCQF(cdbg2, curBlock2, true);
            curBlock2++;
        }
        maxMinimizer =
                // If it's the last block for both of the CQFs
                curBlock1 >= cdbg1.get_numBlocks() and curBlock2 >= cdbg2.get_numBlocks() ? std::max(maxMinimizer1, maxMinimizer2) :
                curBlock1 >= cdbg1.get_numBlocks() ? maxMinimizer2 : // If we've passed the last block for CQF1
                curBlock2 >= cdbg2.get_numBlocks() ? maxMinimizer1 : // If we've passed the last block for CQF2
                std::min(maxMinimizer1, maxMinimizer2);
        console->info("Min minimizer={}, Max minimizer={}", minMinimizer, maxMinimizer);

        //        The output block kmer count should be left for the next set of input blocks
        std::vector<std::pair<uint64_t, uint64_t>> currMinimizerKmers;
        currMinimizerKmers.reserve(cdbg.getCqfSlotCnt());
        for (uint64_t b = minMinimizer; b <= maxMinimizer; b++) {
            // merge the two keys from cqf1 and cqf2
            uint64_t cntr0{0}, cntr1{0}, size0{0}, size1{0};
            if (minimizerKeyList[0][ends[0]-b])
                size0 = minimizerKeyList[0][ends[0]-b]->size();
            if (minimizerKeyList[1][ends[1]-b])
                size1 = minimizerKeyList[1][ends[1]-b]->size();
            if (b % 500 == 0)
                std::cerr << "\rminimizer " << b
                          << " size=" << size0 << "," << size1 << "     ";
            while (cntr0 < size0 and cntr1 < size1) {
                uint64_t key0 = (*minimizerKeyList[0][ends[0]-b])[cntr0];
                uint64_t color0 = (*minimizerColorList[0][ends[0]-b])[cntr0];
                uint64_t key1 = (*minimizerKeyList[1][ends[1]-b])[cntr1];
                uint64_t color1 = (*minimizerColorList[1][ends[1]-b])[cntr1];
                if (key0 < key1) {
                    colorId = get_colorID(std::make_pair(color0, 0));
                    keyObj = KeyObject(key0, 0, colorId);
                    ++cntr0;
                } else if (key0 > key1) {
                    colorId = get_colorID(std::make_pair(0, color1));
                    keyObj = KeyObject(key1, 0, colorId);
                    ++cntr1;
                } else { // it0->first == it1->first
                    colorId = get_colorID(std::make_pair(color0, color1));
                    keyObj = KeyObject(key0, 0, colorId);
                    ++cntr0;
                    ++cntr1;
                }
                currMinimizerKmers.emplace_back(keyObj.key, keyObj.count);
                currMinimizerSlotsCnt++;
                if (colorId > 1) {
                    currMinimizerSlotsCnt += cdbg.getColorIdSlotCnt(colorId) + 1;
                }
                if(colorId <= sampledPairs.size())
                    foundAbundantId++;
            }
            if (cntr0 == size0) {
                while (cntr1 < size1) {
                    uint64_t key1 = (*minimizerKeyList[1][ends[1]-b])[cntr1];
                    uint64_t color1 = (*minimizerColorList[1][ends[1]-b])[cntr1];
                    colorId = get_colorID(std::make_pair(0, color1));
                    currMinimizerKmers.emplace_back(key1, colorId);
                    currMinimizerSlotsCnt++;
                    if (colorId > 1) {
                        currMinimizerSlotsCnt += cdbg.getColorIdSlotCnt(colorId) + 1;
                    }
                    cntr1++;
                    if(colorId <= sampledPairs.size())
                        foundAbundantId++;
                }
            } else {
                while (cntr0 < size0) {
                    uint64_t key0 = (*minimizerKeyList[0][ends[0]-b])[cntr0];
                    uint64_t color0 = (*minimizerColorList[0][ends[0]-b])[cntr0];
                    colorId = get_colorID(std::make_pair(color0, 0));
                    currMinimizerKmers.emplace_back(key0, colorId);
                    currMinimizerSlotsCnt++;
                    if (colorId > 1) {
                        currMinimizerSlotsCnt += cdbg.getColorIdSlotCnt(colorId) + 1;
                    }
                    ++cntr0;
                    if(colorId <= sampledPairs.size())
                        foundAbundantId++;
                }
            }
            for (auto i = 0; i < 2; ++i) {
                minimizerKeyList[i].erase(minimizerKeyList[i].begin()+(ends[i]-b), minimizerKeyList[i].end());
                minimizerColorList[i].erase(minimizerColorList[i].begin()+(ends[i]-b), minimizerColorList[i].end());
                minimizerKeyList[i].shrink_to_fit();
                minimizerColorList[i].shrink_to_fit();
            }
            if (not tmp_kmers.empty() and occupiedSlotsCnt + currMinimizerSlotsCnt > cdbg.getCqfSlotCnt()) {
                std::cerr << "\r";
                console->info("Fill and serialize cqf {} with {} kmers into {} slots up to minimizer {}", outputCQFBlockId, tmp_kmers.size(), occupiedSlotsCnt, b);
                kmerCount+=tmp_kmers.size();
                tbb::parallel_sort(tmp_kmers.begin(), tmp_kmers.end(), [](auto &kv1, auto &kv2) {
                    return kv1.first < kv2.first;
                });
//                std::sort(std::execution::par_unseq, tmp_kmers.begin(), tmp_kmers.end(), [](auto &kv1, auto &kv2) {
//                    return kv1.first < kv2.first;
//                });
                // unique
                tmp_kmers.erase(std::unique(tmp_kmers.begin(), tmp_kmers.end(),
                                           [](auto &kv1, auto &kv2) {
                                               return kv1.first == kv2.first and kv1.second == kv2.second;
                                           }), tmp_kmers.end());
                console->info("Sort-unique done.");
                uint64_t qbits = cdbg.get_qbits();
                console->info("CurrQbits: {}, availableSlotCnt: {}, requiredSlotCnt: {} for {} kmers",
                              qbits, 1ULL << qbits, occupiedSlotsCnt, tmp_kmers.size());
                cdbg.initializeNewCQFBlock(outputCQFBlockId, kbits, qbits, hashmode, seed);
                cdbg.add_kmer2CurDbg(tmp_kmers, 0,  tmp_kmers.size());
                /*int ret;
                do {
                    ret = 0;
                    std::vector<int> rets(numThreads, 0);
                    if (qbits > cdbg.get_qbits()) {
                        console->warn("NOOO!! This thing should not happen very often! Almost never.");
                    }
                    console->info("CurrQbits: {}, availableSlotCnt: {}, requiredSlotCnt: {} for {} kmers",
                                  qbits, 1ULL << qbits, occupiedSlotsCnt, tmp_kmers.size());
                    cdbg.initializeNewCQFBlock(outputCQFBlockId, kbits, qbits, hashmode, seed);
                    std::vector<std::thread> threads;
                    for (uint32_t t = 0; t < numThreads; ++t) {
                        uint64_t s = tmp_kmers.size()*((double)t/numThreads),
                                e = ((double)(t+1)/numThreads)*tmp_kmers.size();
                        threads.emplace_back(std::thread( [&, s, e, t] {
                            std::stringstream ss;
                            ss << "thread " << t << " s=" << s << " e=" << e << " start with reset return value " << rets[t] << "\n";
                            std::cerr << ss.str();
                            rets[t] = cdbg.add_kmer2CurDbg(std::ref(tmp_kmers), s, e);
                        }));
                    }
                    for (uint32_t t = 0; t < numThreads; ++t) {
                        threads[t].join();
                        if (rets[t] == QF_NO_SPACE) {
                            ret = rets[t];
                        }
                    }
                    threads.clear();
                    // increase qbits by 1
                    qbits++;
                } while (ret == QF_NO_SPACE);*/
                cdbg.serializeCurrentCQF();
                console->info("Done filling and storing cqf {}", outputCQFBlockId);
                tmp_kmers.clear();
                tmp_kmers.shrink_to_fit();
//                std::vector<std::remove_reference<decltype(tmp_kmers)>::type::value_type>().swap(tmp_kmers);
                occupiedSlotsCnt = 0;
                outputCQFBlockId++;
            }
            cdbg.minimizerBlock[b] = outputCQFBlockId;
//            std::move(std::execution::par_unseq, currMinimizerKmers.begin(), currMinimizerKmers.end(), tmp_kmers.end());
            tmp_kmers.insert(tmp_kmers.end(), currMinimizerKmers.begin(), currMinimizerKmers.end());
            currMinimizerKmers.clear();
            currMinimizerKmers.shrink_to_fit();
            std::vector<std::remove_reference<decltype(currMinimizerKmers)>::type::value_type>().swap(currMinimizerKmers);
            if (maxMinimizerCnt < currMinimizerSlotsCnt) {
                maxMinimizerCnt = currMinimizerSlotsCnt;
            }
            occupiedSlotsCnt += currMinimizerSlotsCnt;
            currMinimizerSlotsCnt = 0;
        }
        minMinimizer = maxMinimizer + 1;
        std::cerr << "\r";
    }
    for (auto i = 0; i < 2; i++) {
        minimizerKeyList[i].clear();
        minimizerColorList[i].clear();
        minimizerKeyList[i].shrink_to_fit();
        minimizerColorList[i].shrink_to_fit();
    }
    console->info("Maximum number of pairs of key-color per minimizer : {}", maxMinimizerCnt);
    // fill and serialize last cqf block
    if (!tmp_kmers.empty()) {
        console->info("Fill and serialize cqf {} with {} kmers into {} slots as the last cqf block", outputCQFBlockId, tmp_kmers.size(), occupiedSlotsCnt);
        kmerCount+=tmp_kmers.size();
        tbb::parallel_sort(tmp_kmers.begin(), tmp_kmers.end(), [](auto &kv1, auto &kv2) {
            return kv1.first < kv2.first;
        });
//        std::sort(std::execution::par_unseq,tmp_kmers.begin(), tmp_kmers.end(), [](auto &kv1, auto &kv2) {
//            return kv1.first < kv2.first;
//        });
        // unique
        tmp_kmers.erase(std::unique(tmp_kmers.begin(), tmp_kmers.end(),
                                    [](auto &kv1, auto &kv2) {
                                        return kv1.first == kv2.first and kv1.second == kv2.second;
                                    }), tmp_kmers.end());
        console->info("Sort-unique done.");
        auto qbits = static_cast<uint64_t >(ceil(std::log2(occupiedSlotsCnt*100.0/95.0)));
        console->info("Selected qbits for last cqf: {}", qbits);
        console->info("CurrQbits: {}, availableSlotCnt: {}, requiredSlotCnt: {} for {} kmers",
                      qbits, 1ULL << qbits, occupiedSlotsCnt, tmp_kmers.size());
        cdbg.initializeNewCQFBlock(outputCQFBlockId, kbits, qbits, hashmode, seed);
        cdbg.add_kmer2CurDbg(tmp_kmers, 0,  tmp_kmers.size());
/*        int ret;
        do {
            ret = 0;
            std::vector<int> rets(numThreads, 0);
            if (qbits > cdbg.get_qbits()) {
                console->warn("This thing should not happen very often! Almost never. "
                              "CurrQbits: {}, availableSlotCnt: {}, requiredSlotCnt: {} for {} kmers",
                              qbits, 1ULL << qbits, occupiedSlotsCnt, tmp_kmers.size());
            }
            cdbg.initializeNewCQFBlock(outputCQFBlockId, kbits, qbits, hashmode, seed);
            std::vector<std::thread> threads;
            for (uint32_t t = 0; t < numThreads; ++t) {
                uint64_t s = tmp_kmers.size()*((double)t/numThreads),
                e = ((double)(t+1)/numThreads)*tmp_kmers.size();
                threads.emplace_back(std::thread( [&, s, e, t] {
                    std::stringstream ss;
                    ss << "thread " << t << " s=" << s << " e=" << e << " start with reset return value " << rets[t] << "\n";
                    std::cerr << ss.str();
                    rets[t] = cdbg.add_kmer2CurDbg(std::ref(tmp_kmers), s, e);
                }));
            }
            for (uint32_t t = 0; t < numThreads; ++t) {
                threads[t].join();
                if (rets[t] == QF_NO_SPACE) {
                    ret = rets[t];
                }
            }
            threads.clear();
            std::cerr << "Cleared threads\n";
            usleep(10000000);

            // increase qbits by 1
            qbits++;
        } while (ret == QF_NO_SPACE);*/
        console->info("Done constructing last cqf");
        cdbg.serializeCurrentCQF();
        tmp_kmers.clear();
        tmp_kmers.shrink_to_fit();
        occupiedSlotsCnt = 0;
    }


    for (uint64_t m = maxMinimizer+1; m < cdbg.minimizerBlock.size(); m++) {
        cdbg.minimizerBlock[m] = outputCQFBlockId;
    }

    console -> info("Total kmers merged: {}. Time-stamp: {}.", kmerCount, time(nullptr) - start_time_);
    console -> info("Out of {} kmers, {} have color-ids that are sampled earlier.", kmerCount, foundAbundantId);

    tmp_kmers.clear();
    tmp_kmers.shrink_to_fit();
    std::vector<std::remove_reference<decltype(tmp_kmers)>::type::value_type>().swap(tmp_kmers);
    cdbg1.replaceCQFInMemory(invalid);
    cdbg2.replaceCQFInMemory(invalid);
    cdbg.replaceCQFInMemory(invalid);

    auto t_end = time(nullptr);
    console -> info("Merging the CQFs took time {} seconds.", t_end - t_start);
}



template <typename qf_obj, typename key_obj>
uint64_t CQF_merger<qf_obj, key_obj>:: get_colorID(const std::pair<uint64_t, uint64_t> &idPair)
{
    auto it = sampledPairs.find(idPair);
    if(it != sampledPairs.end()) {
        return it->second + 1;
    }
    ColorPair cpair(idPair.first, idPair.second);
    return sampledPairs.size() + colorMph->lookup(cpair) + 1;
}


template <typename qf_obj, typename key_obj>
void CQF_merger<qf_obj, key_obj>::
store_colorID_map()
{
    auto t_start = time(nullptr);

    console -> info("At color-class building (bitvectors concatenation) phase. Time-stamp = {}.",
                    time(nullptr) - start_time_);
    uint64_t writtenPairsCount = 0, cdbg1ColorCnt{cdbg1.get_num_eqclasses()}, cdbg2ColorCnt{cdbg2.get_num_eqclasses()};
    console->info("Writing color pairs and ids into file {}", cdbg.prefix + "newID2oldIDs" );
    std::ofstream output(cdbg.prefix + "newID2oldIDs");
    output.write(reinterpret_cast<char*>(&writtenPairsCount), sizeof(writtenPairsCount));
    output.write(reinterpret_cast<char*>(&cdbg1ColorCnt), sizeof(cdbg1ColorCnt));
    output.write(reinterpret_cast<char*>(&cdbg2ColorCnt), sizeof(cdbg2ColorCnt));
    console->info("# of abundant color IDs: {}", sampledPairs.size());
    // write down the pair and the associated colorID here
    // storing the IDs from 0 (when inserting into CQF, they needed a +1) (0-based)

    for (auto &idpair : sampledPairs) {
        uint64_t colorID = idpair.second;
        auto fs = idpair.first;
        colorIdType c1 = fs.first == 0 ? cdbg1ColorCnt - 1 : fs.first - 1;
        colorIdType c2 = fs.second == 0 ? cdbg2ColorCnt - 1 : fs.second - 1;
        output.write(reinterpret_cast<char*>(&colorID), sizeof(colorID));
        output.write(reinterpret_cast<char*>(&(c1)), sizeof(c1));
        output.write(reinterpret_cast<char*>(&(c2)), sizeof(c2));
    }
    writtenPairsCount = sampledPairs.size();

    // writing down the non-popular color IDs
    std::string colorIdPairFile = cdbg.prefix + EQ_ID_PAIRS_FILE;
    std::ifstream input(colorIdPairFile);
    colorIdType c1, c2;
    while (input >> c1 >> c2) {
        ColorPair cpair(c1, c2);
        uint64_t colorID = sampledPairs.size() + colorMph->lookup(cpair);
        c1 = c1 == 0 ? cdbg1ColorCnt - 1 : c1 - 1;
        c2 = c2 == 0 ? cdbg2ColorCnt - 1: c2 - 1;
        output.write(reinterpret_cast<char*>(&colorID), sizeof(colorID));
        output.write(reinterpret_cast<char*>(&(c1)), sizeof(c1));
        output.write(reinterpret_cast<char*>(&(c2)), sizeof(c2));
        writtenPairsCount++;
    }
    input.close();
    output.seekp(0, std::ios::beg);
    output.write(reinterpret_cast<char*>(&writtenPairsCount), sizeof(writtenPairsCount));
    output.close();

    std::string sysCommand = "rm " + colorIdPairFile;
    system(sysCommand.c_str());
    auto t_end = time(nullptr);
    console -> info("Writing {} color pairs took time {} seconds.", writtenPairsCount, t_end - t_start);
}


template <typename qf_obj, typename key_obj>
void CQF_merger<qf_obj, key_obj>:: serializeRemainingStructures()
{

    store_colorID_map();

    // Serialize the sample-id map.
    std::ofstream outputFile(cdbg.prefix + mantis::SAMPLEID_FILE);

    for(auto idSample : cdbg.sampleid_map)
        outputFile << idSample.first << " " << idSample.second << "\n";

    outputFile.close();

    console -> info("Serialized sample-id mapping.");

    // Serialize minimizer count and associated CQF block
    std::ofstream minfile(cdbg.prefix + mantis::MINIMIZER_FREQ, std::ios::binary);
    minfile.write(reinterpret_cast<char *>(cdbg.minimizerCntr.data()), cdbg.minimizerCntr.size()*sizeof(typename decltype(cdbg.minimizerCntr)::value_type));
    minfile.close();

    minfile.open(cdbg.prefix + mantis::MINIMIZER_BOUNDARY, std::ios::binary);
    minfile.write(reinterpret_cast<char *>(cdbg.minimizerBlock.data()),
                  cdbg.minimizerBlock.size()*sizeof(typename decltype(cdbg.minimizerBlock)::value_type));
    minfile.close();
    console -> info("Serialized minimizer info.");
}

template <typename qf_obj, typename key_obj>
void CQF_merger<qf_obj, key_obj>::merge()
{
    auto t_start = time(nullptr);
    console -> info ("Merging the two CQFs..");
    std::cerr << "\n0. Before starting anything\n";
    auto tillBlock = sample_colorID_pairs(mantis::SAMPLE_SIZE);
    std::cerr << "\n1. After sampling the colorIdPairs and before finding the rest of pairs\n";
    store_colorID_pairs(tillBlock);
    std::cerr << "\n2. After fillDiskBucket before buildMPH\n";
    build_MPHF();
    std::cerr << "\n3. After buildMPH before buildCQF\n";
    build_CQF();
    std::cerr << "\n4. After buildCQF before serialize\n";
    serializeRemainingStructures();
    std::cerr << "\n5. After serialize\n";
    auto t_end = time(nullptr);
    console -> info("CQF merge completed in {} s.", t_end - t_start);
}


template class CQF_merger<SampleObject<CQF<KeyObject> *>, KeyObject>;
