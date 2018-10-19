//
// Created by Fatemeh Almodaresi on 2018-10-18.
//
#include <MantisFS.h>
#include <sstream>
#include "mstQuery.h"
#include "ProgOpts.h"

typedef std::vector<sdsl::rrr_vector < 63>> eqvec;

void loadEqs(spdlog::logger *logger, std::string prefix, eqvec &bvs) {
    std::vector<std::string> eqclass_files =
            mantis::fs::GetFilesExt(prefix.c_str(), mantis::EQCLASS_FILE);

    // sort eqclass_files
    // note to @robP: It terribly statically relies on the format of the input files!!
    std::sort(eqclass_files.begin(), eqclass_files.end(), [logger](std::string &s1, std::string &s2) {
        uint32_t id1, id2;
        std::stringstream ss1(first_part(last_part(s1, '/'), '_'));
        std::stringstream ss2(first_part(last_part(s2, '/'), '_'));
        if ((ss1 >> id1).fail() || !(ss1 >> std::ws).eof() ||
            (ss2 >> id2).fail() || !(ss2 >> std::ws).eof()) {
            logger->error("file name does not start with a number : {}, {}", s1, s2);
        }
        return id1 < id2;
    });
    bvs.reserve(eqclass_files.size());
        uint64_t accumTotalEqCls = 0;
        for (auto &eqfile : eqclass_files) {
            sdsl::rrr_vector<63> bv;
            bvs.push_back(bv);
            sdsl::load_from_file(bvs.back(), eqfile);
        }
    logger->info("Done loading all the equivalence classes.");
}

std::vector<uint64_t> buildColor(eqvec &bvs,
                uint64_t eqid,
                uint64_t num_samples) {
    std::vector<uint64_t> eq;
    eq.reserve(num_samples);
    uint64_t i{0}, bitcnt{0};
    uint64_t idx = eqid / mantis::NUM_BV_BUFFER;
    uint64_t offset = eqid % mantis::NUM_BV_BUFFER;
//std::cerr << eqid << " " << num_samples << " " << idx << " " << offset << "\n";
    while (i<num_samples) {
        bitcnt = std::min(num_samples - i, (uint64_t) 64);
        uint64_t wrd = (bvs[idx]).get_int(offset * num_samples + i, bitcnt);
        for (auto j=0; j < bitcnt; j++) {
            if ((wrd >> j) & 0x01)
                eq.push_back(i+j);
        }
        i += bitcnt;
    }
    return eq;
}


int validate_mst_main(MSTValidateOpts &opt) {
    spdlog::logger *logger = opt.console.get();
    std::string dbg_file(opt.prefix + mantis::CQF_FILE);
    std::string sample_file(opt.prefix + mantis::SAMPLEID_FILE);

    QueryStats queryStats;
    queryStats.numSamples = opt.numSamples;
    logger->info("Number of experiments: {}", queryStats.numSamples);

    logger->info("Loading parentbv, deltabv, and bbv...");
    MSTQuery mstQuery(opt.prefix, opt.k, opt.k, queryStats.numSamples, logger);
    logger->info("Done Loading data structure. Total # of color classes is {}",
                 mstQuery.parentbv.size() - 1);

    logger->info("Loading color classes...");
    eqvec bvs;
    loadEqs(logger, opt.prefix, bvs);
    uint64_t eqCount{0};
    for (auto &bv:bvs) {
        eqCount += bv.size()/opt.numSamples;
    }
    logger->info("Done Loading color classes."
                 "\n\t# of color classes: {}"
                 "\n\t# of Samples: {}", eqCount, opt.numSamples);
    uint64_t cntr{0};
    LRUCacheMap cache_lru(100000);
    for (uint64_t idx = 0; idx < eqCount; idx++) {
        nonstd::optional<uint64_t> dummy{nonstd::nullopt};
        std::vector<uint64_t> newEq = mstQuery.buildColor(idx, queryStats, &cache_lru, nullptr, dummy);
        cache_lru.emplace(idx, newEq);
        std::vector<uint64_t> oldEq = buildColor(bvs, idx, opt.numSamples);
        if (newEq != oldEq) {
            std::cerr << "AAAAA! LOOOSER!!\n";
            std::cerr << "index=" << idx << "\n";
            std::cerr << "new size: " << newEq.size() << " old size: " << oldEq.size() << "\n";
            std::cerr << "new: ";
            for (unsigned long k : newEq) {
                std::cerr << k << " ";
            }
            std::cerr << "\nold: ";
            for (unsigned long k : oldEq) {
                std::cerr << k << " ";
            }
            std::cerr << "\n";
            std::exit(1);
        }
        cntr++;
        if (cntr % 1000000 == 0) {
            std::cerr << "\r" << cntr/1000000 << "M eqs were the same";
        }
    }
    logger->info("\nWOOOOW! Validation passed\n");
}