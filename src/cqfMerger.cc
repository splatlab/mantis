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
    start_time_ = std::time(nullptr);
    kbits = cdbg1.get_current_cqf()->keybits();
    hashmode = cdbg1.get_current_cqf()->hash_mode();
    seed = cdbg1.get_current_cqf()->seed();
    kmerMask = (1ULL << kbits) - 1;

    minimizerKeyColorList[0].resize(cdbg1.minimizerBlock.size());
    minimizerKeyColorList[1].resize(cdbg2.minimizerBlock.size());
    for (auto &m : minimizerKeyColorList[0]) m.reset(new std::vector<std::pair<uint64_t, colorIdType>>());//m->clear();
    for (auto &m : minimizerKeyColorList[1]) m.reset(new std::vector<std::pair<uint64_t, colorIdType>>());//m->clear();
}

template <typename qf_obj, typename key_obj>
uint64_t CQF_merger<qf_obj, key_obj>::
walkBlockedCQF(ColoredDbg<qf_obj, key_obj> &curCdbg, const uint64_t curBlock, bool isSecond) {
    uint64_t minMinimizer{invalid}, maxMinimizer{MAX_MINIMIZER};
    curCdbg.replaceCQFInMemory(curBlock);

    if (curBlock < curCdbg.get_numBlocks()) {
        const CQF<key_obj> *cqf1 = curCdbg.get_current_cqf();
        typename CQF<key_obj>::Iterator it1 = cqf1->begin();
        uint64_t cntr{0}, count{cqf1->dist_elts()};
        auto minmax = curCdbg.getMinMaxMinimizer(curBlock);
        minMinimizer = minmax.first;
        maxMinimizer = minmax.second;

        for (auto i = minMinimizer; i <= maxMinimizer; i++) {
            minimizerKeyColorList[isSecond][i]->reserve(curCdbg.minimizerCntr[i]);
        }
        while (!it1.done()) {
            auto keyval = it1.get_cur_hash();
            auto key = hash_64i(keyval.key, kmerMask);
            auto minimizerPair = curCdbg.findMinimizer(key, kbits);
            if (minimizerPair.first == minimizerPair.second) {
                minimizerPair.second = invalid;
            }
            std::vector<uint64_t> pairs{minimizerPair.first, minimizerPair.second};
            for (auto minimizer : pairs) {
                if (minimizer != invalid and minimizer >= minMinimizer and minimizer <= maxMinimizer) {
                    minimizerKeyColorList[isSecond][minimizer]->emplace_back(keyval.key, keyval.count);
                    cntr++;
                }
            }
            ++it1;
            if (cntr % 10000000 == 0) {
                std::cerr << "\r" << (isSecond?"Second ":"First ") << cntr << " for " << count << " kmers      ";
            }
        }
        std::cerr << "\r" << (isSecond?"Second ":"First ") << cntr << " (main & duplicate) inserts for " << count << " kmers\n";
    }
    return maxMinimizer;
}

template <typename qf_obj, typename key_obj>
uint64_t CQF_merger<qf_obj, key_obj>::
sample_color_id_pairs(uint64_t sampleKmerCount)
{
    auto t_start = time(nullptr);

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

        if (threadCount > 1) {
            std::future<uint64_t> r1 =
                    std::async(std::launch::async, &CQF_merger<qf_obj, key_obj>::walkBlockedCQF, this, std::ref(cdbg1), curBlock, false);
            std::future<uint64_t> r2 =
                    std::async(std::launch::async, &CQF_merger<qf_obj, key_obj>::walkBlockedCQF, this, std::ref(cdbg2), curBlock, true);
            maxMinimizer1 = r1.get();
            maxMinimizer2 = r2.get();
        } else {
            maxMinimizer1 = walkBlockedCQF(cdbg1, curBlock, false);
            maxMinimizer2 = walkBlockedCQF(cdbg2, curBlock, true);
        }
        maxMinimizer =
                // If it's the last block for both of the CQFs
                curBlock >= cdbg1.get_numBlocks() - 1 and curBlock >= cdbg2.get_numBlocks() - 1 ? std::max(maxMinimizer1, maxMinimizer2) :
                curBlock >= cdbg1.get_numBlocks() ? maxMinimizer2 : // If we've passed the last block for CQF1
                curBlock >= cdbg2.get_numBlocks() ? maxMinimizer1 : // If we've passed the last block for CQF2
                std::min(maxMinimizer1, maxMinimizer2);
        std::cerr << "\rMin minimizer=" << minMinimizer << " Max minimizer=" << maxMinimizer << "\n";
        curBlock++;
        for (uint64_t b = minMinimizer; b <= maxMinimizer; b++) {
            // merge the two keys from cqf1 and cqf2
            auto it0 = minimizerKeyColorList[0][b]->begin();
            auto it1 = minimizerKeyColorList[1][b]->begin();
            if (b % 500 == 0)
                std::cerr << "\rminimizer " << b
                          << " size=" << minimizerKeyColorList[0][b]->size() << "," << minimizerKeyColorList[1][b]->size() << "     ";
            while (it0 != minimizerKeyColorList[0][b]->end() and it1 != minimizerKeyColorList[1][b]->end()) {
                if (it0->first < it1->first) {
                    sampledPairs[std::make_pair(it0->second, 0)]++;
                    it0++;
                } else if (it0->first > it1->first) {
                    sampledPairs[std::make_pair(0, it1->second)]++;
                    it1++;
                } else { // it0->first == it1->first
                    sampledPairs[std::make_pair(it0->second, it1->second)]++;
                    it0++;
                    it1++;
                }
                cdbg.minimizerCntr[b]++;
                kmerCount++;
            }
            if (it0 == minimizerKeyColorList[0][b]->end()) {
                while (it1 != minimizerKeyColorList[1][b]->end()) {
                    sampledPairs[std::make_pair(0, it1->second)]++;
                    cdbg.minimizerCntr[b]++;
                    kmerCount++;
                    it1++;
                }
            } else {
                while (it0 != minimizerKeyColorList[0][b]->end()) {
                    sampledPairs[std::make_pair(it0->second, 0)]++;
                    cdbg.minimizerCntr[b]++;
                    kmerCount++;
                    it0++;
                }
            }
            minimizerKeyColorList[0][b].reset(nullptr);
            minimizerKeyColorList[1][b].reset(nullptr);
        }
        minMinimizer = maxMinimizer+1;
        std::cerr << "\r";
    }
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
fill_disk_bucket(uint64_t startingBlock)
{
    auto t_start = time(nullptr);

    console -> info("Writing the non-sampled color-id pairs to disk-files of form ({}). Time-stamp = {}.",
                    EQ_ID_PAIRS_FILE + std::string("_X_Y"), time(nullptr) - start_time_);

    uint64_t writtenPairsCount = 0;


    // force currentCQFBlock to get loaded into memory
    cdbg1.replaceCQFInMemory(invalid);
    cdbg2.replaceCQFInMemory(invalid);

    std::string output = cdbg.prefix + EQ_ID_PAIRS_FILE;
    console -> info("Iterating over the CQFs for the non-sampled color-id pairs starting from block {}"
                    "writing in file \"{}\"", startingBlock, output);
    std::ofstream diskBucket(output, std::ios::out);
    uint64_t curBlock{startingBlock}, kmerCount{0};
    uint64_t maxMinimizer{0}, minMinimizer;
    minMinimizer = curBlock >= cdbg1.get_numBlocks() and curBlock >= cdbg2.get_numBlocks() ? 0 :
            curBlock >= cdbg1.get_numBlocks() ? cdbg2.getMinMaxMinimizer(curBlock).first :
            curBlock >= cdbg2.get_numBlocks() ? cdbg1.getMinMaxMinimizer(curBlock).first :
            std::min(cdbg1.getMinMaxMinimizer(curBlock).first, cdbg2.getMinMaxMinimizer(curBlock).first);

    while(curBlock < cdbg1.get_numBlocks() or curBlock < cdbg2.get_numBlocks()) {
        console->info("Current block={}", curBlock);
        uint64_t maxMinimizer1{0}, maxMinimizer2{0};

        if (threadCount > 1) {
            std::future<uint64_t> r1 = std::async(std::launch::async, &CQF_merger<qf_obj, key_obj>::walkBlockedCQF, this, std::ref(cdbg1), curBlock, false);
            std::future<uint64_t> r2 = std::async(std::launch::async, &CQF_merger<qf_obj, key_obj>::walkBlockedCQF, this, std::ref(cdbg2), curBlock, true);
            maxMinimizer1 = r1.get();
            maxMinimizer2 = r2.get();
        } else {
            maxMinimizer1 = walkBlockedCQF(cdbg1, curBlock, false);
            maxMinimizer2 = walkBlockedCQF(cdbg2, curBlock, true);
        }

        maxMinimizer =
                // If it's the last block for both of the CQFs
                curBlock >= cdbg1.get_numBlocks() - 1 and curBlock >= cdbg2.get_numBlocks() - 1 ? std::max(maxMinimizer1, maxMinimizer2) :
                curBlock >= cdbg1.get_numBlocks() ? maxMinimizer2 : // If we've passed the last block for CQF1
                curBlock >= cdbg2.get_numBlocks() ? maxMinimizer1 : // If we've passed the last block for CQF2
                std::min(maxMinimizer1, maxMinimizer2);
        std::cerr << "\rMin minimizer=" << minMinimizer << " Max minimizer=" << maxMinimizer << "\n";
        curBlock++;
        uint64_t clearedBytes = 0;
        for (uint64_t b = minMinimizer; b <= maxMinimizer; b++) {
            // merge the two keys from cqf1 and cqf2
            auto it0 = minimizerKeyColorList[0][b]->begin();
            auto it1 = minimizerKeyColorList[1][b]->begin();
            if (b % 500 == 0)
                std::cerr << "\rminimizer " << b
                          << " size=" << minimizerKeyColorList[0][b]->size() << "," << minimizerKeyColorList[1][b]->size() << "     ";
            while (it0 != minimizerKeyColorList[0][b]->end() and it1 != minimizerKeyColorList[1][b]->end()) {
                auto colorPair = std::make_pair(it0->second, it1->second);
                if (it0->first < it1->first) {
                    colorPair = std::make_pair(it0->second, 0);
                    it0++;
                } else if (it0->first > it1->first) {
                    colorPair = std::make_pair(0, it1->second);
                    it1++;
                } else { // keys are the same for the two CQFs : it0->first == it1->first
                    it0++;
                    it1++;
                }

                if(sampledPairs.find(colorPair) == sampledPairs.end())
                {
                    diskBucket << colorPair.first << " " << colorPair.second << "\n";
                    writtenPairsCount++;
                }
                cdbg.minimizerCntr[b]++;
                kmerCount++;
            }
            if (it0 == minimizerKeyColorList[0][b]->end()) {
                while (it1 != minimizerKeyColorList[1][b]->end()) {
                    auto colorPair = std::make_pair(0ULL, it1->second);
                    if(sampledPairs.find(colorPair) == sampledPairs.end())
                    {
                        diskBucket << colorPair.first << " " << colorPair.second << "\n";
                        writtenPairsCount++;
                    }
                    it1++;
                    cdbg.minimizerCntr[b]++;
                    kmerCount++;
                }
            } else {
                while (it0 != minimizerKeyColorList[0][b]->end()) {
                    auto colorPair = std::make_pair(it0->second, 0ULL);
                    if(sampledPairs.find(colorPair) == sampledPairs.end())
                    {
                        diskBucket << colorPair.first << " " << colorPair.second << "\n";
                        writtenPairsCount++;
                    }
                    it0++;
                    cdbg.minimizerCntr[b]++;
                    kmerCount++;
                }
            }
            clearedBytes += (minimizerKeyColorList[0][b]->size()+minimizerKeyColorList[1][b]->size());
            minimizerKeyColorList[0][b].reset(nullptr);
            minimizerKeyColorList[1][b].reset(nullptr);
        }
        minMinimizer = maxMinimizer+1;
        std::cerr << "\r";
        std::cerr << "\tcleared MBs in this round: " << (clearedBytes * 12) / std::pow(2, 20) << "\n";
    }
    console -> info("Distinct kmers found {}, color-id pairs written to disk {}. Time-stamp = {}.",
                    kmerCount, writtenPairsCount, time(nullptr) - start_time_);

    diskBucket.flush();
    diskBucket.close();

    auto t_end = time(nullptr);
    console -> info("Filling up the disk-buckets with color-id pairs took time {} seconds.", t_end - t_start);
    filter_disk_buckets();

    return kmerCount;
}


template <typename qf_obj, typename key_obj>
void CQF_merger<qf_obj, key_obj>::
filter_disk_buckets()
{
    if(!system(nullptr))
    {
        console -> error("Command processor does not exist.");
        exit(1);
    }
    uint64_t maxMemoryForSort = 5;

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
build_MPH_tables()
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


    ColorIdPairIterator kb(colorIdPairFile);
    ColorIdPairIterator ke(colorIdPairFile, true);
    auto colorPairIt = boomphf::range(kb, ke);
    colorMph = std::make_unique<boophf_t>(colorIdPairCount, colorPairIt, threadCount, gammaFactor);
    console -> info("Total memory consumed by all the MPH tables = {} MB.", (colorMph->totalBitSize() / 8) / (1024 * 1024));
    auto t_end = time(nullptr);
    console -> info("Building the MPH tables took time {} seconds and {} memory", t_end - t_start, colorMph->totalBitSize());
}

template <typename qf_obj, typename key_obj>
void CQF_merger<qf_obj, key_obj>::build_CQF()
{
    auto t_start = time(nullptr);

    console -> info("At CQFs merging phase. Time-stamp = {}.\n", time(nullptr) - start_time_);

    std::vector<std::pair<uint64_t, uint64_t>> tmp_kmers;
    // reserve 1/3rd more than the threshold as the actual count of kmers is gonna be +-epsilon of the threshold
    tmp_kmers.reserve(block_kmer_threshold + block_kmer_threshold/3);

    uint64_t kmerCount{0}, foundAbundantId{0};
    cdbg1.replaceCQFInMemory(invalid);
    cdbg2.replaceCQFInMemory(invalid);

    uint64_t curBlock{0}, outputCQFBlockId{0}, blockKmerCnt{0};
    for (auto &m :  minimizerKeyColorList[0]) m.reset(new std::vector<std::pair<uint64_t, colorIdType>>());//m.clear();
    for (auto &m :  minimizerKeyColorList[1]) m.reset(new std::vector<std::pair<uint64_t, colorIdType>>());//m.clear();

    uint64_t maxMinimizer{0}, minMinimizer{0};
    uint64_t colorId{0}, epsilon{100};
    KeyObject keyObj;

    cdbg.initializeNewCQFBlock(invalid, kbits, hashmode, seed);
    cdbg.initializeNewCQFBlock(outputCQFBlockId, kbits, hashmode, seed);
    while(curBlock < cdbg1.get_numBlocks() or curBlock < cdbg2.get_numBlocks()) {
        console->info("Current block={}", curBlock);
        uint64_t maxMinimizer1{0}, maxMinimizer2{0};
        if (threadCount > 1) {
            std::future<uint64_t> r1 = std::async(std::launch::async, &CQF_merger<qf_obj, key_obj>::walkBlockedCQF, this, std::ref(cdbg1), curBlock, false);
            std::future<uint64_t> r2 = std::async(std::launch::async, &CQF_merger<qf_obj, key_obj>::walkBlockedCQF, this, std::ref(cdbg2), curBlock, true);
            maxMinimizer1 = r1.get();
            maxMinimizer2 = r2.get();
        } else {
            maxMinimizer1 = walkBlockedCQF(cdbg1, curBlock, false);
            maxMinimizer2 = walkBlockedCQF(cdbg2, curBlock, true);
        }
        maxMinimizer =
                // If it's the last block for both of the CQFs
                curBlock >= cdbg1.get_numBlocks() - 1 and curBlock >= cdbg2.get_numBlocks() - 1 ? std::max(maxMinimizer1, maxMinimizer2) :
                curBlock >= cdbg1.get_numBlocks() ? maxMinimizer2 : // If we've passed the last block for CQF1
                curBlock >= cdbg2.get_numBlocks() ? maxMinimizer1 : // If we've passed the last block for CQF2
                std::min(maxMinimizer1, maxMinimizer2);
        std::cerr << "\rMin minimizer=" << minMinimizer << " Max minimizer=" << maxMinimizer << "\n";
        curBlock++;

        //        The output block kmer count should be left for the next set of input blocks
        for (uint64_t b = minMinimizer; b <= maxMinimizer; b++) {
            // merge the two keys from cqf1 and cqf2
            auto it0 = minimizerKeyColorList[0][b]->begin();
            auto it1 = minimizerKeyColorList[1][b]->begin();
            if (b % 500 == 0)
                std::cerr << "\rminimizer " << b
                          << " size=" << minimizerKeyColorList[0][b]->size() << "," << minimizerKeyColorList[1][b]->size() << "     ";
            while (it0 != minimizerKeyColorList[0][b]->end() and it1 != minimizerKeyColorList[1][b]->end()) {
                if (it0->first < it1->first) {
                    colorId = get_color_id(std::make_pair(it0->second, 0));
                    keyObj = KeyObject(it0->first, 0, colorId);
                    it0++;
                } else if (it0->first > it1->first) {
                    colorId = get_color_id(std::make_pair(0, it1->second));
                    keyObj = KeyObject(it1->first, 0, colorId);
                    it1++;
                } else { // it0->first == it1->first
                    colorId = get_color_id(std::make_pair(it0->second, it1->second));
                    keyObj = KeyObject(it0->first, 0, colorId);
                    it0++;
                    it1++;
                }
                tmp_kmers.emplace_back(keyObj.key, keyObj.count);
                blockKmerCnt++;
                if(colorId <= sampledPairs.size())
                    foundAbundantId++;
            }
            if (it0 == minimizerKeyColorList[0][b]->end()) {
                while (it1 != minimizerKeyColorList[1][b]->end()) {
                    colorId = get_color_id(std::make_pair(0, it1->second));
                    tmp_kmers.emplace_back(it1->first, colorId);
                    blockKmerCnt++;
                    it1++;
                    if(colorId <= sampledPairs.size())
                        foundAbundantId++;
                }
            } else {
                while (it0 != minimizerKeyColorList[0][b]->end()) {
                    colorId = get_color_id(std::make_pair(it0->second, 0));
                    tmp_kmers.emplace_back(it0->first, colorId);
                    blockKmerCnt++;
                    it0++;
                    if(colorId <= sampledPairs.size())
                        foundAbundantId++;
                }
            }
            cdbg.minimizerBlock[b] = outputCQFBlockId;
            minimizerKeyColorList[0][b].reset(nullptr);
            minimizerKeyColorList[1][b].reset(nullptr);
            if (blockKmerCnt and blockKmerCnt > block_kmer_threshold - epsilon) {
                console->info("Fill and serialize cqf {} with {} kmers up to minimizer {}", outputCQFBlockId, blockKmerCnt, b);
                kmerCount += blockKmerCnt;
                blockKmerCnt = 0;
                std::sort(tmp_kmers.begin(), tmp_kmers.end(), [](auto &kv1, auto &kv2){
                    return kv1.first < kv2.first;
                });
                for (auto &kv : tmp_kmers) {
                    keyObj = KeyObject(kv.first, 0, kv.second);
                    cdbg.add_kmer2CurDbg(keyObj, QF_NO_LOCK | QF_KEY_IS_HASH);
                }

                cdbg.serializeCurrentCQF();
                tmp_kmers.clear();

                outputCQFBlockId++;
                cdbg.initializeNewCQFBlock(outputCQFBlockId, kbits, hashmode, seed);
            }
        }
        minMinimizer = maxMinimizer + 1;
        std::cerr << "\r";
    }

    // fill and serialize last cqf block
    if (!tmp_kmers.empty()) {
        console->info("Fill and serialize cqf {} with {} kmers as the last cqf block", outputCQFBlockId, blockKmerCnt);
        kmerCount += blockKmerCnt;
        blockKmerCnt = 0;
        std::sort(tmp_kmers.begin(), tmp_kmers.end(), [](auto &kv1, auto &kv2) {
            return kv1.first < kv2.first;
        });
        for (auto &kv : tmp_kmers) {
            keyObj = KeyObject(kv.first, 0, kv.second);
            cdbg.add_kmer2CurDbg(keyObj, QF_NO_LOCK | QF_KEY_IS_HASH);
        }
        cdbg.serializeCurrentCQF();
        tmp_kmers.clear();
    }


    for (uint64_t m = maxMinimizer+1; m < cdbg.minimizerBlock.size(); m++) {
        cdbg.minimizerBlock[m] = outputCQFBlockId;
    }

    console -> info("Total kmers merged: {}. Time-stamp: {}.", kmerCount, time(nullptr) - start_time_);
    console -> info("Out of {} kmers, {} have color-ids that are sampled earlier.", kmerCount, foundAbundantId);

    cdbg1.replaceCQFInMemory(invalid);
    cdbg2.replaceCQFInMemory(invalid);
    cdbg.replaceCQFInMemory(invalid);

    auto t_end = time(nullptr);
    console -> info("Merging the CQFs took time {} seconds.", t_end - t_start);
}



template <typename qf_obj, typename key_obj>
uint64_t CQF_merger<qf_obj, key_obj>:: get_color_id(const std::pair<uint64_t, uint64_t> &idPair)
{
    auto it = sampledPairs.find(idPair);
    if(it != sampledPairs.end()) {
        return it->second + 1;
    }
    return sampledPairs.size() + colorMph->lookup(idPair) + 1;
}


template <typename qf_obj, typename key_obj>
void CQF_merger<qf_obj, key_obj>::
store_color_pairs()
{
    auto t_start = time(nullptr);

    console -> info("At color-class building (bitvectors concatenation) phase. Time-stamp = {}.",
                    time(nullptr) - start_time_);
    uint64_t writtenPairsCount = 0;
    console->info("Writing color pairs and ids into file {}", cdbg.prefix + "/newID2oldIDs" );
    std::ofstream output(cdbg.prefix + "newID2oldIDs");
    output.write(reinterpret_cast<char*>(&writtenPairsCount), sizeof(writtenPairsCount));

    console->info("# of abundant color IDs: {}", sampledPairs.size());
    // write down the pair and the associated colorID here
    // storing the IDs from 0 (when inserting into CQF, they needed a +1)
    for (auto &idpair : sampledPairs) {
        uint64_t colorID = idpair.second;
        auto fs = idpair.first;
        output.write(reinterpret_cast<char*>(&colorID), sizeof(colorID));
        output.write(reinterpret_cast<char*>(&(fs.first)), sizeof(fs.first));
        output.write(reinterpret_cast<char*>(&(fs.second)), sizeof(fs.second));
    }
    writtenPairsCount = sampledPairs.size();

    // writing down the non-popular color IDs
    std::ifstream input(cdbg.prefix + EQ_ID_PAIRS_FILE);
    std::pair<uint64_t, uint64_t> idPair;
    while (input >> idPair.first >> idPair.second) {
        uint64_t colorID = sampledPairs.size() + colorMph->lookup(idPair);
        output.write(reinterpret_cast<char*>(&colorID), sizeof(colorID));
        output.write(reinterpret_cast<char*>(&(idPair.first)), sizeof(idPair.first));
        output.write(reinterpret_cast<char*>(&(idPair.second)), sizeof(idPair.second));
        writtenPairsCount++;
    }
    input.close();
    output.seekp(0, std::ios::beg);
    output.write(reinterpret_cast<char*>(&writtenPairsCount), sizeof(writtenPairsCount));
    output.close();

    auto t_end = time(nullptr);
    console -> info("Writing {} color pairs took time {} seconds.", writtenPairsCount, t_end - t_start);
}


template <typename qf_obj, typename key_obj>
void CQF_merger<qf_obj, key_obj>:: serializeRemainingStructures()
{

    store_color_pairs();

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
    auto tillBlock = sample_color_id_pairs(mantis::SAMPLE_SIZE);
    fill_disk_bucket(tillBlock);
    build_MPH_tables();
    build_CQF();
    serializeRemainingStructures();
    auto t_end = time(nullptr);
    console -> info("CQF merge completed in {} s.", t_end - t_start);
}


template class CQF_merger<SampleObject<CQF<KeyObject> *>, KeyObject>;