/*
 * =====================================================================================
 *
 *         Author:  Prashant Pandey (), ppandey@cs.stonybrook.edu
 *   Organization:  Stony Brook University
 *
 * =====================================================================================
 */

#include <iostream>
#include <fstream>
#include <utility>
#include <algorithm>
#include <string>
#include <vector>
#include <map>
#include <queue>
#include <set>
#include <unordered_set>
#include <bitset>
#include <cassert>
#include <fstream>

#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdlib.h>
#include <fcntl.h>
#include <unistd.h>
#include <sys/resource.h>
#include <sys/stat.h>
#include <sys/time.h>
#include <sys/mman.h>
#include <openssl/rand.h>

#include "MantisFS.h"
#include "ProgOpts.h"
#include "kmer.h"
#include "coloreddbg.h"
#include "mantisconfig.hpp"
#include "squeakrconfig.h"

#include    <stdlib.h>
#include <mstQuery.h>

/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  main
 *  Description:  
 * =====================================================================================
 */
int
validate_main(ValidateOpts &opt) {

    spdlog::logger *console = opt.console.get();
    float eps = 0.001;
    std::string prefix = opt.prefix;
    if (prefix.back() != '/') {
        prefix += '/';
    }
    // make the output directory if it doesn't exist
    if (!mantis::fs::DirExists(prefix.c_str())) {
        std::cerr << "Couldn't find the input mantis directory.\n";
        std::exit(3);
    }

    // Read the colored dBG
    console->info("Reading colored dbg from disk.");

    ColoredDbg<SampleObject<CQF<KeyObject> *>, KeyObject> cdbg(prefix, MANTIS_DBG_IN_MEMORY);
    cdbg.set_console(console);
    uint64_t kmer_size = cdbg.get_current_cqf()->keybits()/2; // kmer size is in number of bases
    /*console->info("Read colored dbg with {} k-mers and {} color classes",
                  cdbg.get_cqf()->dist_elts(), cdbg.get_num_bitvectors());*/

    std::string query_file = opt.query_file;
    console->info("Reading query kmers from disk.");
    uint64_t total_kmers = 0;
    std::unordered_map<mantis::KmerHash, uint64_t> _dummy_uniqueKmers;
    mantis::QuerySets multi_kmers = Kmer::parse_kmers(query_file.c_str(),
                                                      kmer_size,
                                                      total_kmers,
                                                      false,
                                                      _dummy_uniqueKmers);
    std::vector<uint64_t> queryKmers(multi_kmers.size(), 0);
    uint64_t queryCntr=0;
    for (auto &kmers : multi_kmers) {
        console->info("querySeq {}: {} kmers", queryCntr+1, kmers.size());
        queryKmers[queryCntr] = kmers.size();
        queryCntr++;
    }
    console->info("Total k-mers to query: {}", total_kmers);


    console->info("Query all as a bulk query from MST");
    std::unique_ptr<MSTQuery> mstQuery = std::make_unique<MSTQuery>(prefix, prefix + "querytmps", kmer_size, kmer_size, cdbg.get_num_samples(), console);
    console->info("Done Loading data structure. Total # of color classes is {}",
                  mstQuery->parentbv.size() - 1);
    console->info("Start querying Mantis.");
    LRUCacheMap cache_lru(100000);
    RankScores rs(1);
    QueryStats queryStats;
    queryStats.numSamples = cdbg.get_num_samples();
    uint64_t numOfQueries = mstQuery->parseBulkKmers(query_file, kmer_size);
    console->info("Done reading {} input queries and parsing the kmers.", numOfQueries);
    mstQuery->findSamples(cdbg, cache_lru, &rs, queryStats, numOfQueries);
    mantis::QueryResults result = mstQuery->allQueries; //getResultList(numOfQueries);
    mstQuery.reset(nullptr);
    console->info("Done querying the Mantis index.");



    // Read experiment CQFs
    /*std::ifstream infile(opt.inlist);
    uint64_t num_samples{0};
    if (infile.is_open()) {
        std::string line;
        while (std::getline(infile, line)) { ++num_samples; }
        infile.clear();
        infile.seekg(0, std::ios::beg);
    } else {
        console->error("Input filter list {} does not exist or could not be opened.", opt.inlist);
        std::exit(1);
    }*/

    std::vector<std::pair<uint32_t, std::string>> squeakrFiles;
    // Read experiment CQFs
    std::ifstream sampleFile(prefix + mantis::SAMPLEID_FILE);
    uint64_t num_samples{0};
    if (sampleFile.is_open()) {
        uint64_t order;
        std::string squeakrAdd;
        while (sampleFile >> order >> squeakrAdd) {
            ++num_samples;
            squeakrFiles.emplace_back(order, squeakrAdd);
        }
    } else {
        console->error("Input filter list {} does not exist or could not be opened.", opt.inlist);
        std::exit(1);
    }
    std::sort(squeakrFiles.begin(), squeakrFiles.end(), [](auto &a1, auto &a2) {return a1.first < a2.first;});
    // Query kmers in each experiment CQF in the same order as the Mantis sample file
    uint32_t nqf = 0;
    std::vector<uint64_t> fails;
//    while (infile >> squeakr_file) {
    for (auto & p : squeakrFiles) {
        std::string  squeakr_file = p.second;
        bool fail{false};
        if (!mantis::fs::FileExists(squeakr_file.c_str())) {
            console->error("Squeakr file {} does not exist.", squeakr_file);
            exit(1);
        }
        squeakr::squeakrconfig config;
        int ret = squeakr::read_config(squeakr_file, &config);
        if (ret == squeakr::SQUEAKR_INVALID_VERSION) {
            console->error("Squeakr index version is invalid. Expected: {} Available: {}",
                           squeakr::INDEX_VERSION, config.version);
            exit(1);
        }
        if (kmer_size != config.kmer_size) {
                console->error("Squeakr file {} has a different k-mer size. Expected: {} Available: {}",
                               squeakr_file, kmer_size, config.kmer_size);
                exit(1);
            }
        if (config.cutoff == 1) {
            console->warn("Squeakr file {} is not filtered.", squeakr_file);
        }

        CQF<KeyObject> cqf(squeakr_file, CQF_FREAD);
//        cqf.dump_metadata();
        std::string sample_id = first_part(first_part(last_part(squeakr_file, '/'),
                                                      '.'), '_');
        console->info("Sample {}:{} --> Validating on cdbg", nqf, sample_id);
        queryCntr = 0;
        std::vector<uint64_t> fraction_present(multi_kmers.size(), 0);
        for (auto &kmers : multi_kmers) {
//            console->info("querySeq {}: {} kmers", queryCntr+1, kmers.size());
            for (auto kmer : kmers) {
                KeyObject k(kmer, 0, 0);
                uint64_t count = cqf.query(k, 0);
                if (count > 0) {
                    fraction_present[queryCntr] += 1;
                }
            }
            queryCntr++;
        }

        // Validate the cdbg output
        uint64_t nonEmptyQueries{0}, kmerCnt{0}, mantis_nonEmptyQueries{0}, mantis_kmerCnt{0};
        for (queryCntr = 0; queryCntr < multi_kmers.size(); queryCntr++) {
            if (fraction_present[queryCntr] != result[queryCntr][nqf]) {
                console->error("Failed for query {}, original CQF {}, cdbg {} out of {}",
                              queryCntr, fraction_present[queryCntr], result[queryCntr][nqf], queryKmers[queryCntr]);
                fail = true;
                //abort();
            }
            if (fraction_present[queryCntr] > 0) {
                nonEmptyQueries += 1;
                kmerCnt += fraction_present[queryCntr];
            }
            if (result[queryCntr][nqf] > 0) {
                mantis_nonEmptyQueries += 1;
                mantis_kmerCnt += result[queryCntr][nqf];
            }
        }
        /*if (!cqfs.front().check_similarity(&cqfs.back())) {
            console->error("Passed Squeakr files are not similar.", squeakr_file);
            exit(1);
        }*/
        if (fail) {
            console->error("Failed --> Sample {}: {}! Avg kmer count: Mantis: {}, Squeakrs: {}", nqf, sample_id,
                           static_cast<uint64_t >(mantis_kmerCnt/(mantis_nonEmptyQueries+eps)),
                           static_cast<uint64_t >(kmerCnt/(nonEmptyQueries+eps)));
            fails.push_back(nqf);
//            std::exit(3);
        } else {
            console->info("Passed --> {} non-zero intersected queries with avg found kmer count={} ", nonEmptyQueries,
                          static_cast<uint64_t >(kmerCnt/(nonEmptyQueries+eps)));
        }
        nqf++;
    }

    console->info("");
    if (not fails.empty()) {
        console->error("Mantis validation failed!");
        console->info("Sample IDs that Mantis validation failed on:");
        for (auto f : fails) {
            std::cerr << f << " ";
        }
        std::cerr << "\n";
        console->error("Mantis validation failed!");
    }
    console->info("YaaaY");
    console->info("Mantis validation passed!");

#if 0
    // This is x-axis
    // For one query set
    // Divide aggregate by the total kmers in the query to get the fraction.
    uint64_t cnt = 0;
    std::unordered_map<uint64_t, float> fraction_present = ground_truth[0];
    for (uint64_t i = 0; i < nqf; i++)
        fraction_present[i] = (fraction_present[i] / (float)multi_kmers[cnt++].size()) * 100;

    std::vector<std::vector<std::string>> buckets;
    for (auto it : fraction_present)
        buckets[it.second / 5].push_back(cdbg.get_sample(it.first));
#endif

    return EXIT_SUCCESS;
}                /* ----------  end of function main  ---------- */

