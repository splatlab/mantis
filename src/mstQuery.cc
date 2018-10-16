//
// Created by Fatemeh Almodaresi on 2018-10-15.
//
#include <fstream>
#include <vector>
#include <CLI/Timer.hpp>

#include "ProgOpts.h"
#include "kmer.h"
#include "mstQuery.h"

void MSTQuery::loadIdx(std::string indexDir) {
    sdsl::load_from_file(parentbv, indexDir + mantis::PARENTBV_FILE);
    sdsl::load_from_file(deltabv, indexDir + mantis::DELTABV_FILE);
    sdsl::load_from_file(bbv, indexDir + mantis::BOUNDARYBV_FILE);
    sbbv = sdsl::bit_vector::select_1_type(&bbv);
    zero = parentbv.size() - 1; // maximum color id which
    logger->info("Loaded the new color class index");
    logger->info("\t--> parent size: ", parentbv.size());
    logger->info("\t--> delta size: ", deltabv.size());
    logger->info("\t--> boundary size: ", bbv.size());
}

std::vector<uint64_t> MSTQuery::buildColor(uint64_t eqid, QueryStats &queryStats,
                                           LRUCacheMap *lru_cache,
                                           RankScores *rs,
                                           nonstd::optional<uint64_t> &toDecode, // output param.  Also decode these
                                           bool all = true) {
    (void) rs;
    std::vector<uint32_t> flips(numSamples);
    std::vector<uint32_t> xorflips(numSamples, 0);
    uint64_t i{eqid}, from{0}, to{0};
    int64_t height{0};
    auto &froms = queryStats.buffer;
    froms.clear();
    queryStats.totEqcls++;
    bool foundCache = false;
    uint32_t iparent = parentbv[i];
    while (iparent != i) {
        if (lru_cache and lru_cache->contains(i)) {
            const auto &vs = (*lru_cache)[i];
            for (auto v : vs) {
                xorflips[v] = 1;
            }
            queryStats.cacheCntr++;
            foundCache = true;
            break;
        }
        from = (i > 0) ? (sbbv(i) + 1) : 0;
        froms.push_back(from);

        if (queryStats.trySample) {
            auto &occ = queryStats.numOcc[iparent];
            ++occ;
            if ((!toDecode) and
                (occ > 10) and
                (height > 10) and
                (lru_cache and
                 !lru_cache->contains(iparent))) {
                toDecode = iparent;
            }
        }
        i = iparent;
        iparent = parentbv[i];
        ++queryStats.totSel;
        ++height;
    }
    if (!foundCache and i != zero) {
        from = (i > 0) ? (sbbv(i) + 1) : 0;
        froms.push_back(from);
        ++queryStats.totSel;
        queryStats.rootedNonZero++;
        ++height;
    }
    uint64_t pctr{0};
    for (auto f : froms) {
        bool found = false;
        uint64_t wrd{0};
        uint64_t offset{0};
        auto start = f;
        do {
            wrd = bbv.get_int(start, 64);
            for (uint64_t j = 0; j < 64; j++) {
                flips[deltabv[start + j]] ^= 0x01;
                if ((wrd >> j) & 0x01) {
                    found = true;
                    break;
                }
            }
            start += 64;
        } while (!found);
    }

    if (!all) { // return the indices of set bits
        std::vector<uint64_t> eq;
        eq.reserve(numWrds);
        uint64_t one = 1;
        for (i = 0; i < numSamples; i++) {
            if (flips[i] ^ xorflips[i]) {
                eq.push_back(i);
            }
        }
        return eq;
    }
    std::vector<uint64_t> eq(numWrds);
    uint64_t one = 1;
    for (i = 0; i < numSamples; i++) {
        if (flips[i] ^ xorflips[i]) {
            uint64_t idx = i / 64;
            eq[idx] = eq[idx] | (one << (i % 64));
        }
    }
    return eq;
}

mantis::QueryResult MSTQuery::findSamples(const mantis::QuerySet &kmers,
                                          CQF<KeyObject> &dbg,
                                          MSTQuery &mstQuery,
                                          LRUCacheMap &lru_cache,
                                          RankScores *rs,
                                          QueryStats &queryStats) {
    std::unordered_map<uint64_t, uint64_t> query_eqclass_map;
    for (auto k : kmers) {
        KeyObject key(k, 0, 0);
        uint64_t eqclass = dbg.query(key, 0);
        if (eqclass)
            query_eqclass_map[eqclass] += 1;
    }

    mantis::QueryResult sample_map(queryStats.numSamples, 0);
    size_t numPerLevel = 10;
    nonstd::optional<uint64_t> toDecode{nonstd::nullopt};
    nonstd::optional<uint64_t> dummy{nonstd::nullopt};

    for (auto &it : query_eqclass_map) {
        auto eqclass_id = it.first - 1;
        auto count = it.second;

        std::vector<uint64_t> setbits;
        if (lru_cache.contains(eqclass_id)) {
            setbits = lru_cache[eqclass_id];//.get(eqclass_id);
            queryStats.cacheCntr++;
        } else {
            queryStats.noCacheCntr++;
            toDecode.reset();
            dummy.reset();
            queryStats.trySample = (queryStats.noCacheCntr % 10 == 0);
            setbits = buildColor(eqclass_id, queryStats, &lru_cache, rs, toDecode, false);
            lru_cache.emplace(eqclass_id, setbits);
            if ((queryStats.trySample) and toDecode) {
                auto s = buildColor(*toDecode, queryStats, nullptr, nullptr, dummy, false);
                lru_cache.emplace(*toDecode, s);
            }
        }
        for (auto sb : setbits) {
            sample_map[sb] += count;
        }

        ++queryStats.globalQueryNum;
    }

    return sample_map;
}


QueryStats output_results(mantis::QuerySets &multi_kmers,
                          CQF<KeyObject> &dbg,
                          MSTQuery &mstQuery,
                          std::ofstream &opfile,
                          std::vector<std::string> &sampleNames,
                          LRUCacheMap &cache_lru) {
    mantis::QueryResults qres;
    QueryStats queryStats;
    for (auto &kmers : multi_kmers) {
        CLI::AutoTimer timer{"Query time ", CLI::Timer::Big};
        opfile << "seq" << queryStats.cnt++ << '\t' << kmers.size() << '\n';
        mantis::QueryResult result = mstQuery.findSamples(kmers,
                                                          dbg, mstQuery, cache_lru,
                                                          nullptr, queryStats);
        for (uint64_t i = 0; i < result.size(); i++) {
            if (result[i] > 0)
                opfile << sampleNames[i] << '\t' << result[i] << '\n';
        }
    }
    return queryStats;
}

QueryStats output_results_json(mantis::QuerySets &multi_kmers,
                               CQF<KeyObject> &dbg,
                               MSTQuery &mstQuery,
                               std::ofstream &opfile,
                               std::vector<std::string> &sampleNames,
                               LRUCacheMap &cache_lru) {
    mantis::QueryResults qres;
    uint32_t cnt = 0;
    CLI::AutoTimer timer{"Query time ", CLI::Timer::Big};
    opfile << "[\n";
    size_t qctr{0};
    size_t nquery{multi_kmers.size()};
    QueryStats queryStats;
    for (auto &kmers : multi_kmers) {
        CLI::AutoTimer timer{"Query time ", CLI::Timer::Big};
        opfile << "{ \"qnum\": " << queryStats.cnt++ << ",  \"num_kmers\": "
               << kmers.size() << ", \"res\": {\n";
        mantis::QueryResult result = mstQuery.findSamples(kmers,
                                                          dbg, mstQuery, cache_lru,
                                                          nullptr, queryStats);
        uint64_t kmerCntr = 0;
        for (auto it = result.begin(); it != result.end(); ++it) {
            if (*it > 0)
                opfile << " \"" << sampleNames[kmerCntr] << "\": " << *it;
            if (std::next(it) != result.end()) {
                opfile << ",\n";
            }
            kmerCntr++;
        }
        opfile << "}}";
        if (qctr < nquery - 1) {
            opfile << ",";
        }
        opfile << "\n";
        ++qctr;
    }
    opfile << "]\n";

    return queryStats;
}

std::vector<std::string> loadSampleFile(const std::string &sampleFileAddr) {
    std::vector<std::string> sampleNames;
    std::ifstream sampleFile(sampleFileAddr);
    uint64_t id;
    std::string name;
    while (sampleFile >> id >> name) {
        if (id >= sampleNames.size()) {
            sampleNames.resize(id + 1);
        }
        sampleNames[id] = name;
    }
    return sampleNames;
}

/*
 * ===  FUNCTION  =============================================================
 *         Name:  main
 *  Description:
 * ============================================================================
 */
int mst_query_main(QueryOpts &opt) {
    QueryStats queryStats;

    spdlog::logger *logger = opt.console.get();
    std::string dbg_file(opt.prefix + mantis::CQF_FILE);
    std::string sample_file(opt.prefix + mantis::SAMPLEID_FILE);

    std::vector<std::string> sampleNames = loadSampleFile(sample_file);
    queryStats.numSamples = sampleNames.size();
    logger->info("Number of experiments: {}", queryStats.numSamples);

    logger->info("Loading cqf...");
    CQF<KeyObject> cqf(dbg_file, CQF_FREAD);
    auto k = cqf.keybits() / 2;
    logger->info("Done loading cqf. k is {}", k);

    logger->info("Loading color classes...");
    MSTQuery mstQuery(opt.prefix, queryStats.numSamples, logger);
    logger->info("Done Loading color classes. Total # of color classes is {}",
                 mstQuery.parentbv.size() - 1);

    logger->info("Loading query file...");
    uint64_t total_kmers = 0;
    mantis::QuerySets multi_kmers = Kmer::parse_kmers(opt.query_file.c_str(),
                                                      k,
                                                      total_kmers);
    logger->info("Done loading query file : # of seqs: ", multi_kmers.size());


    logger->info("Querying the colored dbg.");
    std::ofstream opfile(opt.output);
    LRUCacheMap cache_lru(100000);
    RankScores rs(1);
    if (opt.use_json) {
        queryStats = output_results_json(multi_kmers, cqf, mstQuery, opfile, sampleNames, cache_lru);
    } else {
        queryStats = output_results(multi_kmers, cqf, mstQuery, opfile, sampleNames, cache_lru);
    }
    opfile.close();
    logger->info("Writing done.");

    logger->info("cache was used {} times and not used {} times",
                 queryStats.cacheCntr, queryStats.noCacheCntr);
    logger->info("total selects = {}, time per select = {}",
                 queryStats.totSel, queryStats.selectTime.count() / queryStats.totSel);
    logger->info("total # of queries = {}, total # of queries rooted at a non-zero node = {}",
                 queryStats.totEqcls, queryStats.rootedNonZero);
/*
    logger->info("select time was {}s, flip time was {}",
            queryStats.selectTime.count(), queryStats.flipTime.count());
*/
/*for (auto &kv : queryStats.numOcc) {
    std::cout << kv.first << '\t' << kv.second << '\n';
}*/
    return EXIT_SUCCESS;
}