//
// Created by Fatemeh Almodaresi on 2020-03-12.
//

#include "cdBG_merger.h"



template<typename qf_obj, typename key_obj>
CdBG_merger<qf_obj, key_obj>::
CdBG_merger(ColoredDbg<qf_obj, key_obj> &&cdbg1in, ColoredDbg<qf_obj, key_obj> &&cdbg2in,
            ColoredDbg<qf_obj, key_obj> &&cdbgOut)
{
    start_time_ = std::time(nullptr);

    if(cdbg1.get_eq_class_file_count() >= cdbg2.get_eq_class_file_count())
        this -> cdbg1 = std::move(cdbg1in), this -> cdbg2 = std::move(cdbg2in);
    else
    {
        this -> cdbg1 = std::move(cdbg2in), this -> cdbg2 = std::move(cdbg1in);

//		console -> info("Mantis indices are swapped.");
    }

    cdbg = std::move(cdbgOut);
    numCCPerBuffer1 = mantis::BV_BUF_LEN / cdbg1.num_samples;
    numCCPerBuffer2 = mantis::BV_BUF_LEN / cdbg2.num_samples;
    numCCPerBuffer = mantis::BV_BUF_LEN / (cdbg1.num_samples + cdbg2.num_samples);
    kbits = cdbg1.get_current_cqf()->keybits();
    hashmode = cdbg1.get_current_cqf()->hash_mode();
    seed = cdbg1.get_current_cqf()->seed();
    kmerMask = (1ULL << kbits) - 1;
}

template <typename qf_obj, typename key_obj>
uint64_t CdBG_merger<qf_obj, key_obj>::
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
uint64_t CdBG_merger<qf_obj, key_obj>::
sample_color_id_pairs(uint64_t sampleKmerCount)
{
    auto t_start = time(nullptr);

    console -> info("Sampling at least {} kmers. Time-stamp = {}.",
                    sampleKmerCount, time(nullptr) - start_time_);

    // force currentCQFBlock to get loaded into memory
//	cdbg1.replaceCQFInMemory(invalid);
//	cdbg2.replaceCQFInMemory(invalid);

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
                    std::async(std::launch::async, &CdBG_merger<qf_obj, key_obj>::walkBlockedCQF, this, std::ref(cdbg1), curBlock, false);
            std::future<uint64_t> r2 =
                    std::async(std::launch::async, &CdBG_merger<qf_obj, key_obj>::walkBlockedCQF, this, std::ref(cdbg2), curBlock, true);
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

    uint64_t colorId = 1;
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
void CdBG_merger<qf_obj, key_obj>::
init_disk_buckets()
{
    const uint64_t fileCount1 = cdbg1.get_eq_class_file_count(),
            fileCount2 = cdbg2.get_eq_class_file_count();

    console -> info("Initializing {} x {} disk-buckets.", fileCount1 + 1, fileCount2 + 1);

    diskBucket.resize(fileCount1 + 1),
            bucketSize.resize(fileCount1 + 1),
            cumulativeBucketSize.resize(fileCount1 + 1),
            MPH.resize(fileCount1 + 1);

    for(uint64_t i = 0; i <= fileCount1; ++i)
    {
        diskBucket[i].resize(fileCount2 + 1),
                bucketSize[i].resize(fileCount2 + 1),
                cumulativeBucketSize[i].resize(fileCount2 + 1),
                MPH[i].resize(fileCount2 + 1);

        for(uint64_t j = 0; j <= fileCount2; ++j)
            diskBucket[i][j] = std::ofstream(cdbg.prefix + TEMP_DIR + EQ_ID_PAIRS_FILE + "_" +
                                             std::to_string(i) + "_" + std::to_string(j)),
            bucketSize[i][j] = 0,
            cumulativeBucketSize[i][j] = 0,
            MPH[i][j] = NULL;
    }

    console -> info("{} x {} disk-buckets  and associated data structures initialized.",
                    fileCount1 + 1, fileCount2 + 1);
}

template <typename qf_obj, typename key_obj>
uint64_t CdBG_merger<qf_obj, key_obj>::
fill_disk_buckets(uint64_t startingBlock)
{
    auto t_start = time(nullptr);

    console -> info("Writing the non-sampled color-id pairs to disk-files of form ({}). Time-stamp = {}.",
                    TEMP_DIR + EQ_ID_PAIRS_FILE + std::string("_X_Y"), time(nullptr) - start_time_);

    uint64_t writtenPairsCount = 0;

    console -> info("Iterating over the CQFs for the non-sampled color-id pairs starting from block {}", startingBlock);

    // force currentCQFBlock to get loaded into memory
    cdbg1.replaceCQFInMemory(invalid);
    cdbg2.replaceCQFInMemory(invalid);

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
            std::future<uint64_t> r1 = std::async(std::launch::async, &CdBG_merger<qf_obj, key_obj>::walkBlockedCQF, this, std::ref(cdbg1), curBlock, false);
            std::future<uint64_t> r2 = std::async(std::launch::async, &CdBG_merger<qf_obj, key_obj>::walkBlockedCQF, this, std::ref(cdbg2), curBlock, true);
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
                colorCount1 = std::max(colorCount1, static_cast<uint64_t >(it0->second));
                colorCount2 = std::max(colorCount2, static_cast<uint64_t >(it1->second));
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
                    add_color_id_pair(colorPair.first, colorPair.second, diskBucket);
                    writtenPairsCount++;
                }
                cdbg.minimizerCntr[b]++;
                kmerCount++;
            }
            if (it0 == minimizerKeyColorList[0][b]->end()) {
                while (it1 != minimizerKeyColorList[1][b]->end()) {
                    colorCount2 = std::max(colorCount2, static_cast<uint64_t >(it1->second));
                    auto colorPair = std::make_pair(0ULL, it1->second);
                    if(sampledPairs.find(colorPair) == sampledPairs.end())
                    {
                        add_color_id_pair(colorPair.first, colorPair.second, diskBucket);
                        writtenPairsCount++;
                    }
                    it1++;
                    cdbg.minimizerCntr[b]++;
                    kmerCount++;
                }
            } else {
                while (it0 != minimizerKeyColorList[0][b]->end()) {
                    colorCount1 = std::max(colorCount1, static_cast<uint64_t >(it0->second));
                    auto colorPair = std::make_pair(it0->second, 0ULL);
                    if(sampledPairs.find(colorPair) == sampledPairs.end())
                    {
                        add_color_id_pair(colorPair.first, colorPair.second, diskBucket);
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


    // Flush and close the disk-bucket files.
    for(int i = 0; i <= cdbg1.get_eq_class_file_count(); ++i)
        for(int j = 0; j <= cdbg2.get_eq_class_file_count(); ++j)
        {
            diskBucket[i][j].flush();
            diskBucket[i][j].close();
        }


    auto t_end = time(nullptr);
    console -> info("Filling up the disk-buckets with color-id pairs took time {} seconds.", t_end - t_start);


    return kmerCount;
}

template <typename qf_obj, typename key_obj>
void CdBG_merger<qf_obj, key_obj>::
add_color_id_pair(const uint64_t colorID1, const uint64_t colorID2,
                  std::vector<std::vector<std::ofstream>> &diskBucket)
{
    // TODO: Add faster file-write mechanism.

    const uint64_t row = (colorID1 ? (colorID1 - 1) / numCCPerBuffer1 + 1 : 0),//mantis::NUM_BV_BUFFER + 1 : 0),
            col = (colorID2 ? (colorID2 - 1) / numCCPerBuffer2 + 1 : 0);//mantis::NUM_BV_BUFFER + 1 : 0);

    diskBucket[row][col] << colorID1 << " " << colorID2 << "\n";
//    diskBucket[0][0] << colorID1 << " " << colorID2 << "\n";
}

template <typename qf_obj, typename key_obj>
uint64_t CdBG_merger<qf_obj, key_obj>::
filter_disk_buckets()
{
    if(!system(NULL))
    {
        console -> error("Command processor does not exist.");
        exit(1);
    }


    auto t_start = time(nullptr);

    console -> info("Filtering out the unique eq-id pairs from files {} with {} threads. Time-stamp = {}",
                    TEMP_DIR + EQ_ID_PAIRS_FILE + "_X_Y", threadCount, time(nullptr) - start_time_);


    uint64_t colorClassCount = 0;
    const uint64_t fileCount1 = cdbg1.get_eq_class_file_count(),
            fileCount2 = cdbg2.get_eq_class_file_count();

    uint maxMemoryForSort = 5;

    for(int i = 0; i <= fileCount1; ++i)
        for(int j = 0; j <= fileCount2; ++j)
        {
            std::string diskBucket = cdbg.prefix + TEMP_DIR + EQ_ID_PAIRS_FILE +
                                     "_" + std::to_string(i) + "_" + std::to_string(j);

            std::string sysCommand = "sort -u";
            sysCommand += " --parallel=" + std::to_string(threadCount);
            sysCommand += " -S " + std::to_string(maxMemoryForSort) + "G";
            sysCommand += " -o " + diskBucket + " " + diskBucket;

            console -> info("System command used:\n{}", sysCommand);

            system(sysCommand.c_str());


            std::string lineCountFile = cdbg.prefix + TEMP_DIR + ID_PAIR_COUNT_FILE;
            sysCommand = "wc -l " + diskBucket + " | egrep -o \"[0-9]+ \" > " + lineCountFile;
            system(sysCommand.c_str());

            std::ifstream lineCount(lineCountFile);

            lineCount >> bucketSize[i][j];
            lineCount.close();

            colorClassCount += bucketSize[i][j];

            console -> info("Filtered {} unique color-id pairs at disk-bucket {}. Time-stamp = {}.",
                            bucketSize[i][j], diskBucket, time(nullptr) - start_time_);
        }


    console -> info("Count of unique color-id pairs = {}. Time-stamp = {}.",
                    colorClassCount, time(nullptr) - start_time_);


    auto t_end = time(nullptr);
    console -> info("Filtering the unique color-id pairs took time {} seconds.", t_end - t_start);

    return colorClassCount;
}

template <typename qf_obj, typename key_obj>
void CdBG_merger<qf_obj, key_obj>::
build_MPH_tables()
{
    auto t_start = time(nullptr);

    console -> info("Building an MPH (Minimal Perfect Hash) table per disk-bucket. Time-stamp = {}.",
                    time(nullptr) - start_time_);


    const uint64_t fileCount1 = cdbg1.get_eq_class_file_count(),
            fileCount2 = cdbg2.get_eq_class_file_count();
    uint64_t mphMemory = 0;


    for(uint64_t i = 0; i <= fileCount1; ++i)
        for(uint64_t j = 0; j <= fileCount2; ++j)
            if(!bucketSize[i][j])
                console -> info("Bucket ({}, {}) is empty. No MPH table is built.", i, j);
            else
            {
                std::vector<std::pair<uint64_t, uint64_t>> colorIdPair;
                colorIdPair.reserve(bucketSize[i][j]);


                console -> info("Loading the unique color-id pairs from bucket ({}, {}) into memory. Time-stamp = {}.",
                                i, j, time(nullptr) - start_time_);

                std::ifstream input(cdbg.prefix + TEMP_DIR + EQ_ID_PAIRS_FILE +
                                    "_" + std::to_string(i) + "_" + std::to_string(j));
                uint64_t id1, id2;

                // TODO: Add faster file-read mechanism.
                while(input >> id1 >> id2)
                    colorIdPair.emplace_back(id1, id2);

                input.close();

                console -> info("Loaded {} color-id pairs into memory. Time-stamp = {}.", colorIdPair.size(),
                                time(nullptr) - start_time_);


                console -> info("Constructing a BooPHF with {} elements using {} threads. Time-stamp = {}.",
                                colorIdPair.size(), threadCount, time(nullptr) - start_time_);

                // Lowest bit/elem is achieved with gamma=1, higher values lead to larger mphf but faster construction/query.
                double gammaFactor = 2.0;	// gamma = 2 is a good tradeoff (leads to approx 3.7 bits/key).

                // Build the mphf.
                MPH[i][j] = new boophf_t(colorIdPair.size(), colorIdPair, threadCount, gammaFactor);
                mphMemory += MPH[i][j] -> totalBitSize();


                console -> info("For bucket ({}, {}) BooPHF constructed perfect hash for {} keys; total memory = {} MB, bits/elem : {}. Time-stamp = {}.",
                                i, j, colorIdPair.size(), (MPH[i][j] -> totalBitSize() / 8) / (1024 * 1024),
                                (double)(MPH[i][j] -> totalBitSize()) / colorIdPair.size(),
                                time(nullptr) - start_time_);


                colorIdPair.clear();
                colorIdPair.shrink_to_fit();
            }

    console -> info("Total memory consumed by all the MPH tables = {} MB.", (mphMemory / 8) / (1024 * 1024));


    auto t_end = time(nullptr);
    console -> info("Building the MPH tables took time {} seconds.", t_end - t_start);
}

template <typename qf_obj, typename key_obj>
void CdBG_merger<qf_obj, key_obj>::build_CQF()
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

    std::cerr << "\n\n\n0SLEEEP before initializing and loading the output cqf";
    usleep(10000000);
    cdbg.initializeNewCQFBlock(invalid, kbits, hashmode, seed);
    cdbg.initializeNewCQFBlock(outputCQFBlockId, kbits, hashmode, seed);
    while(curBlock < cdbg1.get_numBlocks() or curBlock < cdbg2.get_numBlocks()) {
        console->info("Current block={}", curBlock);
        uint64_t maxMinimizer1{0}, maxMinimizer2{0};
        if (threadCount > 1) {
            std::future<uint64_t> r1 = std::async(std::launch::async, &CdBG_merger<qf_obj, key_obj>::walkBlockedCQF, this, std::ref(cdbg1), curBlock, false);
            std::future<uint64_t> r2 = std::async(std::launch::async, &CdBG_merger<qf_obj, key_obj>::walkBlockedCQF, this, std::ref(cdbg2), curBlock, true);
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

        std::cerr << "\n\n\n2SLEEEP loading the two and filling minimizers\n";
        usleep(10000000);

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
                console->info("3SLEEEP before sorting");
                usleep(10000000);

                std::sort(tmp_kmers.begin(), tmp_kmers.end(), [](auto &kv1, auto &kv2){
                    return kv1.first < kv2.first;
                });
                for (auto &kv : tmp_kmers) {
                    keyObj = KeyObject(kv.first, 0, kv.second);
                    cdbg.add_kmer2CurDbg(keyObj, QF_NO_LOCK | QF_KEY_IS_HASH);
                }

                cdbg.serializeCurrentCQF();
                tmp_kmers.clear();
                console->info("4SLEEEP after sorting, before initializing a new one");
                usleep(10000000);

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


    auto t_end = time(nullptr);
    console -> info("Merging the CQFs took time {} seconds.", t_end - t_start);

    console->info("5SLEEEP befor freeing the memory");
    usleep(10000000);

    cdbg1.replaceCQFInMemory(invalid);
    cdbg2.replaceCQFInMemory(invalid);
    cdbg.replaceCQFInMemory(invalid);
    console->info("6SLEEEP after freeing the memory. Before the start of merging MSTs");
    usleep(10000000);

}



template <typename qf_obj, typename key_obj>
uint64_t CdBG_merger<qf_obj, key_obj>:: get_color_id(const std::pair<uint64_t, uint64_t> &idPair)
{
    auto it = sampledPairs.find(idPair);
    if(it != sampledPairs.end()) {
        return it->second;
    }

    const uint64_t row = (idPair.first ? (idPair.first - 1) / numCCPerBuffer1 + 1 : 0),//mantis::NUM_BV_BUFFER + 1 : 0),
            col = (idPair.second ? (idPair.second - 1) / numCCPerBuffer2 + 1 : 0);//mantis::NUM_BV_BUFFER + 1 : 0);
    if (row >= MPH.size() or col >= MPH[row].size()) {
        std::cerr << "Shouldn't happen! row or column greater than MPH size\n";
        std::cerr << row << " " << col << " " << idPair.first << " " << idPair.second << "\n";
    }
    return cumulativeBucketSize[row][col] + MPH[row][col]->lookup(idPair) + 1;
}


template <typename qf_obj, typename key_obj>
void CdBG_merger<qf_obj, key_obj>:: serialize_cqf_and_sampleid_list()
{
    // Serialize the CQF.
    cdbg.serializeCurrentCQF();
    console -> info("Serialized CQF.");


    // Serialize the sample-id map.
    std::ofstream outputFile(cdbg.prefix + mantis::SAMPLEID_FILE);

    for(auto idSample : cdbg.sampleid_map)
        outputFile << idSample.first << " " << idSample.second << "\n";

    outputFile.close();

    console -> info("Serialized sample-id mapping.");

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
void CdBG_merger<qf_obj, key_obj>::
calc_mst_stats(ColoredDbg<qf_obj, key_obj> &cdbg1, ColoredDbg<qf_obj, key_obj> &cdbg2, std::string& dir1, std::string& dir2)
{
    auto t_start = time(nullptr);

    console -> info("At color-class building (bitvectors concatenation) phase. Time-stamp = {}.",
                    time(nullptr) - start_time_);


    const uint64_t bucketCnt1 = cdbg1.get_eq_class_file_count(), bucketCnt2 = cdbg2.get_eq_class_file_count();

    uint64_t writtenPairsCount = 0;
    BitVectorRRR bitVec1, bitVec2;

    std::cerr << dir1 + mantis::PARENTBV_FILE << "\n";
    std::cerr << dir2 + mantis::PARENTBV_FILE << "\n";
    sdsl::int_vector<> parentbv1;
    sdsl::load_from_file(parentbv1, dir1 + "/" + mantis::PARENTBV_FILE);
    sdsl::int_vector<> parentbv2;
    sdsl::load_from_file(parentbv2, dir2 + "/" + mantis::PARENTBV_FILE);

    std::cerr << "loaded the files\n";
    spp::sparse_hash_map<std::pair<uint64_t, uint64_t >, bool, Custom_Pair_Hasher> isVirtual;

    std::cerr << "parents sizes: " << parentbv1.size() << " " << parentbv2.size() << "\n";
    uint64_t zero1 = parentbv1.size() - 1;
    uint64_t zero2 = parentbv2.size() - 1;

    std::cerr << "# of buckets: " << bucketCnt1 << " " << bucketCnt2 << "\n";
    for (uint64_t i = 0; i <= bucketCnt1; ++i)
    {
        std::cerr << "\n" << i;
        for(uint64_t j = 0; j <= bucketCnt2; ++j) {
            std::cerr << "\n" << j << "\n";
            std::ifstream input(cdbg.prefix + TEMP_DIR + EQ_ID_PAIRS_FILE +
                                "_" + std::to_string(i) + "_" + std::to_string(j));
            std::pair<uint64_t, uint64_t> idPair;
            uint64_t cntr{0};
            while (input >> idPair.first >> idPair.second) {
//				std::cerr << "reading the file: " << idPair.first << " " << idPair.second << "\n";
                if (!idPair.first) idPair.first = zero1;
                else idPair.first--;
                if (!idPair.second) idPair.second = zero2;
                else idPair.second--;
                // doesn't matter if the idPair is in the map or not, we want to set the value to false
                isVirtual[idPair] = false;
                uint64_t pid1{idPair.first}, cid1{idPair.first}, pid2{idPair.second}, cid2{idPair.second};
                std::cerr << "\r" << cntr++;
                do {
                    cid1 = pid1;
                    cid2 = pid2;
                    pid1 = parentbv1[cid1];
                    pid2 = parentbv2[cid2];
                    auto ppair = std::make_pair(pid1, pid2);
                    if (isVirtual.find(ppair) == isVirtual.end()) {
                        isVirtual[ppair] = true;
                    }
                } while (pid1 != cid1 and pid2 != cid2);
            }
        }
    }

    uint64_t actualCnt{0}, virtualCnt{0};
    for (auto& kv : isVirtual) {
        if (kv.second) virtualCnt++;
        else actualCnt++;
    }

    std::cerr << "actual cnt: " << actualCnt << " virtual cnt: " << virtualCnt << "\n";
    std::exit(1);
}

template <typename qf_obj, typename key_obj>
uint64_t CdBG_merger<qf_obj, key_obj>::
store_abundant_color_pairs(std::ofstream& output)
{
    console->info("# of abundant color IDs: {}", sampledPairs.size());
    for (auto &idpair : sampledPairs) {
        uint64_t colorID = idpair.second-1;
        auto fs = idpair.first;
        output.write(reinterpret_cast<char*>(&colorID), sizeof(colorID));
        output.write(reinterpret_cast<char*>(&(fs.first)), sizeof(fs.first));
        output.write(reinterpret_cast<char*>(&(fs.second)), sizeof(fs.second));
//		std::cerr << colorID << " " << fs.first << " " << fs.second << "\n";
    }
    return sampledPairs.size();
}


template <typename qf_obj, typename key_obj>
void CdBG_merger<qf_obj, key_obj>::
store_color_pairs(ColoredDbg<qf_obj, key_obj> &cdbg1, ColoredDbg<qf_obj, key_obj> &cdbg2,
                  uint64_t& numColorBuffers)
{
    auto t_start = time(nullptr);

    console -> info("At color-class building (bitvectors concatenation) phase. Time-stamp = {}.",
                    time(nullptr) - start_time_);


    const uint64_t fileCount1 = cdbg1.get_eq_class_file_count(), fileCount2 = cdbg2.get_eq_class_file_count();

    console->info("EQClass Count for cdbg1 and 2: {}, {}", fileCount1, fileCount2);
    uint64_t writtenPairsCount = 0;
    //BitVectorRRR bitVec1, bitVec2;

    console->info("Writing color pairs and ids into file {}", cdbg.prefix + "/newID2oldIDs" );
    std::ofstream output(cdbg.prefix + "newID2oldIDs");
    output.write(reinterpret_cast<char*>(&writtenPairsCount), sizeof(writtenPairsCount));
    writtenPairsCount = store_abundant_color_pairs(output);

    for(uint64_t i = 0; i <= fileCount1; ++i) {

        for (uint64_t j = 0; j <= fileCount2; ++j) {

            console->info("At bucket ({}, {}), size = {}. Time-stamp = {}.", i, j, bucketSize[i][j],
                          time(nullptr) - start_time_);

            std::ifstream input(cdbg.prefix + TEMP_DIR + EQ_ID_PAIRS_FILE +
                                "_" + std::to_string(i) + "_" + std::to_string(j));
            std::pair<uint64_t, uint64_t> idPair;
            cumulativeBucketSize[i][j] = writtenPairsCount;
            while (input >> idPair.first >> idPair.second) {
                /// @fatemeh write down the pair and the associated colorID here
                uint64_t colorID = cumulativeBucketSize[i][j] + MPH[i][j]->lookup(idPair);// + 1;
                output.write(reinterpret_cast<char*>(&colorID), sizeof(colorID));
                output.write(reinterpret_cast<char*>(&(idPair.first)), sizeof(idPair.first));
                output.write(reinterpret_cast<char*>(&(idPair.second)), sizeof(idPair.second));
//				std::cerr << colorID << " " << idPair.first << " " << idPair.second << "\n";

                writtenPairsCount++;
            }
            input.close();
        }
    }
    output.seekp(0, std::ios::beg);
    output.write(reinterpret_cast<char*>(&writtenPairsCount), sizeof(writtenPairsCount));
    output.close();
    numColorBuffers = writtenPairsCount/numCCPerBuffer + 1;//mantis::NUM_BV_BUFFER + 1;

    auto t_end = time(nullptr);
    console -> info("Writing {} color pairs took time {} seconds.", writtenPairsCount, t_end - t_start);
}



template <typename qf_obj, typename key_obj>
void CdBG_merger<qf_obj, key_obj>::merge()
{

    auto t_start = time(nullptr);
    console -> info ("Splitting output minimizers into blocks based on sum of the two input minimizers");
    minimizerKeyColorList[0].resize(cdbg1.minimizerBlock.size());
    minimizerKeyColorList[1].resize(cdbg2.minimizerBlock.size());
    for (auto &m : minimizerKeyColorList[0]) m.reset(new std::vector<std::pair<uint64_t, colorIdType>>());//m->clear();
    for (auto &m : minimizerKeyColorList[1]) m.reset(new std::vector<std::pair<uint64_t, colorIdType>>());//m->clear();
    console -> info ("Done dividing. Time-stamp = {}.\n", time(nullptr) - t_start);

    t_start = time(nullptr);
    console -> info ("Merge starting. Time-stamp = {}.\n", time(nullptr) - start_time_);


    // Make the temporary directory if it doesn't exist.
    std::string tempDir = cdbg.prefix + TEMP_DIR;

    if(!mantis::fs::DirExists(tempDir.c_str()))
        mantis::fs::MakeDir(tempDir.c_str());

    // Check to see if the temporary directory exists now.
    if(!mantis::fs::DirExists(tempDir.c_str()))
    {
        console -> error("Temporary directory {} could not be created.", tempDir);
        exit(1);
    }


    auto tillBlock = sample_color_id_pairs(mantis::SAMPLE_SIZE);

    init_disk_buckets();
    fill_disk_buckets(tillBlock);
    uint64_t colorClassCount = sampledPairs.size() + filter_disk_buckets();
    build_MPH_tables();
    uint64_t num_colorBuffers = 1;
    store_color_pairs(cdbg1, cdbg2, num_colorBuffers);
    console->info("# of color buffers is {}", num_colorBuffers);

    build_CQF();

    // Remove the temporary directory.
    std::string sysCommand = "rm -rf " + tempDir;
    console -> info("Removing the temporary directory. System command used:\n{}", sysCommand);

    system(sysCommand.c_str());


    auto t_end = time(nullptr);

    console -> info("CQF merge completed in {} s.", t_end - t_start);

    auto t_mst_start = time(nullptr);
//    uint64_t num_colorBuffers = 1;
    console->info("{}, {}", cdbg1.prefix, cdbg2.prefix);
    MSTMerger mst(cdbg.prefix, console, threadCount, cdbg1.prefix, cdbg2.prefix, num_colorBuffers);
    console->info("MST Initiated. Now merging the two MSTs..");
    mst.mergeMSTs();
    t_end = time(nullptr);
    console->info("MST merge completed in {} s.", t_end - t_mst_start);
    serialize_cqf_and_sampleid_list();
    t_end = time(nullptr);
    console->info("Total merge time is {} s", t_end - t_start);
}


template class CdBG_merger<SampleObject<CQF<KeyObject> *>, KeyObject>;