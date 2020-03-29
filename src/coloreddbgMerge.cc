/*
 * ============================================================================
 *
 *         Author:  Jamshed Khan (jamshed@cs.umd.edu)
 *   Organization:  University of Maryland, College Park
 *
 * ============================================================================
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

#include "sparsepp/spp.h"
#include "tsl/sparse_map.h"

#include "MantisFS.h"
#include "ProgOpts.h"
#include "coloreddbg.h"
#include "cqfMerger.h"
#include "squeakrconfig.h"
#include "mantisconfig.hpp"


bool data_exists(std::string &dir, spdlog::logger *console) {
    if (mantis::fs::GetFilesExt(dir.c_str(), mantis::CQF_FILE).empty()) {
        console->error("CQF file {} does not exist in input directory {}.", mantis::CQF_FILE, dir);
        return false;
    }

    if (!mantis::fs::FileExists((dir + mantis::SAMPLEID_FILE).c_str())) {
        console->error("Sample-ID list file {} does not exist in input directory {}.",
                       mantis::SAMPLEID_FILE, dir);
        return false;
    }

    /*if(mantis::fs::GetFilesExt(dir.c_str(), mantis::EQCLASS_FILE).empty())
    {
        console -> error("No equivalence-class file with extension {} exists in input directory {}.",
                        mantis::EQCLASS_FILE, dir);
        return false;
    }*/


    return true;
}


int merge_main(MergeOpts &opt) {
    spdlog::logger *console = opt.console.get();

    if (opt.dir1.back() != '/')    // Make sure it is a full directory.
        opt.dir1 += '/';
    // Make sure if the first input directory exists.
    if (!mantis::fs::DirExists(opt.dir1.c_str())) {
        console->error("Input directory {} does not exist.", opt.dir1);
        exit(1);
    }

    if (opt.dir2.back() != '/')    // Make sure it is a full directory.
        opt.dir2 += '/';
    // Make sure if the second input directory exists.
    if (!mantis::fs::DirExists(opt.dir2.c_str())) {
        console->error("Input directory {} does not exist.", opt.dir2);
        exit(1);
    }

    // Check if all the required data exist in the input directories.
    if (!data_exists(opt.dir1, console) || !data_exists(opt.dir2, console))
        exit(1);

    auto t_start = time(nullptr);
    console->info("Merge starting ...");

    // merging two CQFs
    console->info("Merging the two CQFs...");
    auto cqfMerger = std::make_unique<CQF_merger<SampleObject<CQF<KeyObject> *>, KeyObject>>(opt.dir1, opt.dir2,
                                                                                             opt.out, console,
                                                                                             opt.threadCount);
	cqfMerger->merge();
    cqfMerger.reset(nullptr); // make the memory back to almost 0 to start the next section
    // merging two MSTs
    console->info("Merging the two MSTs...");
//    usleep(30000000);
    MSTMerger mst(opt.out, console, opt.threadCount, opt.dir1, opt.dir2);
    mst.mergeMSTs();

    auto t_end = time(nullptr);
    console->info("Total merge time is {} s", t_end - t_start);
    return EXIT_SUCCESS;
}


int validate_merge_main(ValidateMergeOpts &opt) {
    spdlog::logger *console = opt.console.get();

    console->info("Checking correctness of the merged colored dBG.");

    std::string correctResPrefix(opt.correctRes);
    if (correctResPrefix.back() != '/')
        correctResPrefix += '/';

    std::string mergedResPrefix(opt.mergeRes);
    if (mergedResPrefix.back() != '/')
        mergedResPrefix += '/';


    std::string corrCQFfile(correctResPrefix + mantis::CQF_FILE);
    std::vector<std::string> corrColorClassFiles = mantis::fs::GetFilesExt(correctResPrefix.c_str(),
                                                                           mantis::EQCLASS_FILE);
    std::string corrSampleListFile(correctResPrefix + mantis::SAMPLEID_FILE);


    console->info("Loading the correct CdBG into memory.");
    ColoredDbg<SampleObject<CQF<KeyObject> *>, KeyObject> corrCdBG(corrCQFfile, corrColorClassFiles,
                                                                   corrSampleListFile, MANTIS_DBG_IN_MEMORY);

    console->info("Loaded the correct CdBG; it has {} k-mers, {} color-class files, and {} color-classes.",
                  corrCdBG.get_cqf()->dist_elts(), corrColorClassFiles.size(),
                  corrCdBG.get_num_bitvectors());


    std::string mergedCQFfile = mergedResPrefix + mantis::CQF_FILE;
    std::vector<std::string> mergedColorClassFiles = mantis::fs::GetFilesExt(mergedResPrefix.c_str(),
                                                                             mantis::EQCLASS_FILE);
    std::string mergedSampleListFile = mergedResPrefix + mantis::SAMPLEID_FILE;

    console->info("Loading the merged CdBG into memory.");
    ColoredDbg<SampleObject<CQF<KeyObject> *>, KeyObject> mergedCdBG(mergedCQFfile, mergedColorClassFiles,
                                                                     mergedSampleListFile, MANTIS_DBG_IN_MEMORY);

    console->info("Loaded the merged CdBG; it has {} k-mers, {} color-class files, and {} color-classes.",
                  mergedCdBG.get_cqf()->dist_elts(), mergedColorClassFiles.size(),
                  mergedCdBG.get_num_bitvectors());


    if (corrCdBG.get_cqf()->dist_elts() != mergedCdBG.get_cqf()->dist_elts() ||
        corrColorClassFiles.size() != mergedColorClassFiles.size() ||
        corrCdBG.get_num_bitvectors() != mergedCdBG.get_num_bitvectors() ||
        corrCdBG.get_num_samples() != mergedCdBG.get_num_samples()) {
        console->error("Mismatching meta-info.");
        exit(1);
    } else
        console->info("Meta information matches.");


    auto eqclass_m = mergedCdBG.get_eqclasses(), eqclass_c = corrCdBG.get_eqclasses();

    // Linear scan over the CQFs
    auto it_m = mergedCdBG.get_cqf()->begin(), it_c = corrCdBG.get_cqf()->begin();
    uint64_t kmerCount = 0;
    const uint64_t PROGRESS_STEP = 10000000;


    while (!it_m.done() && !it_c.done()) {
        KeyObject entry_m = it_m.get_cur_hash(), entry_c = it_c.get_cur_hash();

        if (entry_m.key != entry_c.key) {
            console->error("Mismatching k-mers found.");
            exit(1);
        }


        uint64_t eqID_m = entry_m.count, eqID_c = entry_c.count;
        std::unordered_set<std::string> sampleSet_m, sampleSet_c;
        const uint64_t wordLen = 64;


        uint64_t sampleCount = mergedCdBG.get_num_samples();
        uint64_t numCCPerBuffer = mantis::BV_BUF_LEN / sampleCount;

        // Fetch the set of samples in the merged CdBG containing the k-mer of this iteration (entry_m.key).
        uint64_t bucketIdx_m = (eqID_m - 1) / numCCPerBuffer;
        uint64_t offset_m = ((eqID_m - 1) % numCCPerBuffer) * sampleCount;

        for (uint32_t wordCount = 0; wordCount <= sampleCount / wordLen; ++wordCount) {
            uint64_t readLen = std::min(wordLen, sampleCount - wordCount * wordLen);
            uint64_t word = eqclass_m[bucketIdx_m].get_int(offset_m, readLen);

            for (uint32_t bitIdx = 0, sampleID = wordCount * wordLen; bitIdx < readLen; ++bitIdx, ++sampleID)
                if ((word >> bitIdx) & 0x01)
                    sampleSet_m.insert(mergedCdBG.get_sample(sampleID));

            offset_m += readLen;
        }


        // Fetch the set of samples in the correct CdBG containing the k-mer of this iteration (entry_c.key).
        uint64_t bucketIdx_c = (eqID_c - 1) / numCCPerBuffer;//mantis::NUM_BV_BUFFER;
        uint64_t offset_c = ((eqID_c - 1) % numCCPerBuffer) * sampleCount;

        for (uint32_t wordCount = 0; wordCount <= sampleCount / wordLen; ++wordCount) {
            uint64_t readLen = std::min(wordLen, sampleCount - wordCount * wordLen);
            uint64_t word = eqclass_c[bucketIdx_c].get_int(offset_c, readLen);

            for (uint32_t bitIdx = 0, sampleID = wordCount * wordLen; bitIdx < readLen; ++bitIdx, ++sampleID)
                if ((word >> bitIdx) & 0x01)
                    sampleSet_c.insert(corrCdBG.get_sample(sampleID));

            offset_c += readLen;
        }

        if (sampleSet_c.size() != sampleSet_m.size()) {
            console->error("For one or more k-mers, the sample set sizes are found to be different between the CdBGs.");
            exit(1);
        }


        for (auto sample : sampleSet_c)
            if (sampleSet_m.find(sample) == sampleSet_m.end()) {
                console->error("Sample sets mismatch for one or more kmers.");
                exit(1);
            }


        ++it_m, ++it_c, kmerCount++;

        if (kmerCount % PROGRESS_STEP == 0)
            console->info("{}M k-mers matched.", kmerCount * 10 / PROGRESS_STEP);
    }


    console->info("CQF sizes = {}, matching k-mers found {}.", corrCdBG.get_cqf()->dist_elts(), kmerCount);
    console->info("Merged CdBG has correct CQF, and correct color-classes for all k-mers.");

    return EXIT_SUCCESS;
}