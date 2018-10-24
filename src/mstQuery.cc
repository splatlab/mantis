//
// Created by Fatemeh Almodaresi on 2018-10-15.
//
#include <fstream>
#include <vector>
#include <CLI/Timer.hpp>
#include <canonicalKmer.h>

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
    logger->info("\t--> parent size: {}", parentbv.size());
    logger->info("\t--> delta size: {}", deltabv.size());
    logger->info("\t--> boundary size: {}", bbv.size());
}

std::vector<uint64_t> MSTQuery::buildColor(uint64_t eqid, QueryStats &queryStats,
                                           LRUCacheMap *lru_cache,
                                           RankScores *rs,
                                           nonstd::optional<uint64_t> &toDecode) {
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
        //std::cerr << i << " " << iparent << "\n";
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
    //std::cerr << "\n\nnext step:\n";
    for (auto f : froms) {
        bool found = false;
        uint64_t wrd{0};
        uint64_t offset{0};
        auto start = f;
        //std::cerr << "\n" << start << ": ";
        do {
            wrd = bbv.get_int(start, 64);
            for (uint64_t j = 0; j < 64; j++) {
                //std::cerr << deltabv[start + j] << " ";
                flips[deltabv[start + j]] ^= 0x01;
                if ((wrd >> j) & 0x01) {
                    found = true;
                    break;
                }
            }
            start += 64;
        } while (!found);
    }

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

void MSTQuery::findSamples(CQF<KeyObject> &dbg,
                           LRUCacheMap &lru_cache,
                           RankScores *rs,
                           QueryStats &queryStats) {
    mantis::EqMap query_eqclass_map;
    std::unordered_set<uint64_t> query_eqclass_set;
//    std::cerr << "\n\nkmer2cidMap size: " << kmer2cidMap.size() << "\n\n";
    for (auto &kv : kmer2cidMap) {
        KeyObject key(kv.first, 0, 0);
//        std::cerr << std::string(dna::canonical_kmer(dbg.keybits()/2, key.key)) << "\n";
        uint64_t eqclass = dbg.query(key, 0);
        if (eqclass) {
            kv.second = eqclass - 1;
            query_eqclass_set.insert(eqclass - 1);
        }
    }

    mantis::QueryResult sample_map(queryStats.numSamples, 0);
    nonstd::optional<uint64_t> toDecode{nonstd::nullopt};
    nonstd::optional<uint64_t> dummy{nonstd::nullopt};

    for (auto &it : query_eqclass_set) {
        uint64_t eqclass_id = it;

        std::vector<uint64_t> setbits;
        if (lru_cache.contains(eqclass_id)) {
            setbits = lru_cache[eqclass_id];//.get(eqclass_id);
            queryStats.cacheCntr++;
        } else {
            queryStats.noCacheCntr++;
            toDecode.reset();
            dummy.reset();
            queryStats.trySample = (queryStats.noCacheCntr % 10 == 0);
            setbits = buildColor(eqclass_id, queryStats, &lru_cache, rs, toDecode);
            lru_cache.emplace(eqclass_id, setbits);
            if ((queryStats.trySample) and toDecode) {
                auto s = buildColor(*toDecode, queryStats, nullptr, nullptr, dummy);
                lru_cache.emplace(*toDecode, s);
            }
        }
        cid2expMap[eqclass_id] = setbits;
        /*for (auto sb : setbits) {
            sample_map[sb] += count;
        }

        ++queryStats.globalQueryNum;*/
    }
}


uint64_t MSTQuery::parseKmers(std::string filename, uint64_t kmer_size) {
    std::ifstream ipfile(filename);
    std::string read;
    uint64_t numOfQueries{0};
    CLI::AutoTimer timer{"First round going over the file ", CLI::Timer::Big};
    while (ipfile >> read) {
        numOfQueries++;
        bool done = false;
        while (!done and read.length() >= kmer_size) {
            uint64_t first = 0;
            uint64_t first_rev = 0;
            uint64_t item = 0;
            bool allNuclValid = true;
            for (uint32_t i = 0; i < kmer_size; i++) { //First kmer
                uint8_t curr = Kmer::map_base(read[i]);
                if (curr > DNA_MAP::G) { // 'N' is encountered
                    if (i + 1 < read.length())
                        read = read.substr(i + 1);
                    else
                        done = true; // you've reached the end of a sequence (last nucl. is N)
                    allNuclValid = false;
                    break;
                }
                first = first | curr;
                first = first << 2;
            }
            if (allNuclValid) {
                first = first >> 2;
                first_rev = static_cast<uint64_t >(Kmer::reverse_complement(first, kmer_size));

                //cout << "kmer: "; cout << int_to_str(first);
                //cout << " reverse-comp: "; cout << int_to_str(first_rev) << endl;

                if (Kmer::compare_kmers(first, first_rev))
                    item = first;
                else
                    item = first_rev;
                kmer2cidMap[item] = 0;

                uint64_t next = (first << 2) & BITMASK(2 * kmer_size);
                uint64_t next_rev = first_rev >> 2;
                uint64_t i = 0;
                for (i = kmer_size; i < read.length(); i++) { //next kmers
                    //cout << "K: " << read.substr(i-K+1,K) << endl;
                    uint8_t curr = Kmer::map_base(read[i]);
                    if (curr > DNA_MAP::G) { // 'N' is encountered
                        if (i + 1 < read.length())
                            read = read.substr(i + 1);
                        else
                            done = true;
                        break;
                    }
                    next |= curr;
                    auto tmp = static_cast<uint64_t>(Kmer::reverse_complement_base(curr));
                    tmp <<= (kmer_size * 2 - 2);
                    next_rev = next_rev | tmp;
                    if (Kmer::compare_kmers(next, next_rev))
                        item = next;
                    else
                        item = next_rev;
                    kmer2cidMap[item] = 0;
                    next = (next << 2) & BITMASK(2 * kmer_size);
                    next_rev = next_rev >> 2;
                }
                if (i == read.length()) done = true;
            }
        }
    }
    return numOfQueries;
}

mantis::QueryResult MSTQuery::convertIndexK2QueryK(std::string &read) {
    mantis::QueryResult res(numSamples, 0);
    auto requiredCnt = static_cast<uint8_t>(queryK - indexK + 1);
    auto queryIndxKDiff = queryK - indexK + 1;
//    std::cerr << "required count: " << requiredCnt << " num samples: " << numSamples << "\n";
    bool done = false;
//    std::cerr << "\n" << read << ": ";
    while (!done and read.length() >= queryK) {
//        std::cerr << read << "\n";
        uint64_t idx2replace = 0;
        std::vector<uint8_t> samples(numSamples);
        std::vector<std::vector<bool>> pastKmers(queryIndxKDiff);
        for (auto &v : pastKmers) {
            v.resize(numSamples, false);
        }
        uint64_t first = 0;
        uint64_t first_rev = 0;
        uint64_t item = 0;
        bool allNuclValid = true;
        for (uint32_t i = 0; i < indexK; i++) { //First kmer
            uint8_t curr = Kmer::map_base(read[i]);
            if (curr > DNA_MAP::G) { // 'N' is encountered
                if (i + 1 < read.length())
                    read = read.substr(i + 1);
                else
                    done = true; // you've reached the end of a sequence (last nucl. is N)
                allNuclValid = false;
                break;
            }
            first = first | curr;
            first = first << 2;
        }
        if (allNuclValid) {
            first = first >> 2;
            first_rev = static_cast<uint64_t >(Kmer::reverse_complement(first, indexK));
            item = Kmer::compare_kmers(first, first_rev) ? first : first_rev;
            pastKmers[idx2replace].assign(pastKmers[idx2replace].size(), false);
//            std::cerr << cid2expMap[kmer2cidMap[item]].size() << ", ";
            for (auto &c : cid2expMap[kmer2cidMap[item]]) {
//                std::cerr << c << " " << pastKmers[idx2replace][c] << " ";
                samples[c]++;
                pastKmers[idx2replace][c] = true;
//                std::cerr << pastKmers[idx2replace][c] << ",";
            }
//            std::cerr << "\n";

            uint64_t next = (first << 2) & BITMASK(2 * indexK);
            uint64_t next_rev = first_rev >> 2;

            for (uint64_t i = indexK; i < queryK; i++) { //next kmers
                idx2replace = ++idx2replace % queryIndxKDiff;
                uint8_t curr = Kmer::map_base(read[i]);
                if (curr > DNA_MAP::G) { // 'N' is encountered
                    if (i + 1 < read.length())
                        read = read.substr(i + 1);
                    else
                        done = true;
                    allNuclValid = false;
                    break;
                }
                next |= curr;
                auto tmp = static_cast<uint64_t>(Kmer::reverse_complement_base(curr));
                tmp <<= (indexK * 2 - 2);
                next_rev = next_rev | tmp;
                item = Kmer::compare_kmers(next, next_rev) ? next : next_rev;
//                std::cerr << cid2expMap[kmer2cidMap[item]].size() << ", ";
                for (auto &c : cid2expMap[kmer2cidMap[item]]) {
                    samples[c]++;
                    pastKmers[idx2replace][c] = true;
                }
                next = (next << 2) & BITMASK(2 * indexK);
                next_rev = next_rev >> 2;
            }
            if (allNuclValid) {
                idx2replace = 0;
                uint64_t i = 0;
                for (i = queryK; i < read.length(); i++) { //next kmers
                    for (auto c = 0; c < samples.size(); c++) {
                        if (samples[c] == requiredCnt) {
                            res[c]++;
                        }
                        if (pastKmers[idx2replace][c]) {
                            samples[c]--;
                        }
                        pastKmers[idx2replace][c] = false;
                    }
//                    std::cerr << i << "\n";
                    uint8_t curr = Kmer::map_base(read[i]);
                    if (curr > DNA_MAP::G) { // 'N' is encountered
                        if (i + 1 < read.length())
                            read = read.substr(i + 1);
                        else
                            done = true;
                        break;
                    }
                    next |= curr;
                    auto tmp = static_cast<uint64_t>(Kmer::reverse_complement_base(curr));
                    tmp <<= (indexK * 2 - 2);
                    next_rev = next_rev | tmp;
                    item = Kmer::compare_kmers(next, next_rev) ? next : next_rev;
//                    std::cerr << cid2expMap[kmer2cidMap[item]].size() << ", ";
                    for (auto &c : cid2expMap[kmer2cidMap[item]]) {
                        samples[c]++;
                        pastKmers[idx2replace][c] = true;
                    }
                    idx2replace = ++idx2replace % queryIndxKDiff;
                    next = (next << 2) & BITMASK(2 * indexK);
                    next_rev = next_rev >> 2;
                }
//                std::cerr << "next\n";
                if (i == read.length())  {
                    for (auto c = 0; c < samples.size(); c++) {
                        if (samples[c] == requiredCnt) {
                            res[c]++;
                        }
                    }
                    done = true;
                }
            }
        }
    }
    return res;
}

void output_results(std::string &queryFile,
                    MSTQuery &mstQuery,
                    std::ofstream &opfile,
                    std::vector<std::string> &sampleNames,
                    QueryStats &queryStats) {
    std::ifstream ipfile(queryFile);
    std::string read;
    CLI::AutoTimer timer{"Second round going over the file + query time ", CLI::Timer::Big};
    while (ipfile >> read) {
        opfile << "seq" << queryStats.cnt++ << '\t' << read.length() << '\n';
        mantis::QueryResult result = mstQuery.convertIndexK2QueryK(read);
//        std::cerr << "result.size(): " << result.size() << "\n";
        for (uint64_t i = 0; i < result.size(); i++) {
            if (result[i] > 0) {
//                std::cerr << i << "\n";
                opfile << sampleNames[i] << '\t' << result[i] << '\n';
            }
        }
//        std::cerr << "next read\n";
    }
}

void output_results_json(std::string &queryFile,
                         MSTQuery &mstQuery,
                         std::ofstream &opfile,
                         std::vector<std::string> &sampleNames,
                         QueryStats &queryStats,
                         uint64_t nquery) {
    std::ifstream ipfile(queryFile);
    std::string read;
    opfile << "[\n";
    uint64_t qctr{0};
    CLI::AutoTimer timer{"Query time ", CLI::Timer::Big};
    while (ipfile >> read) {
        opfile << "{ \"qnum\": " << queryStats.cnt++ << ",  \"num_kmers\": "
               << read.length() << ", \"res\": {\n";
        mantis::QueryResult result = mstQuery.convertIndexK2QueryK(read);
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
    uint32_t queryK = opt.k;
    QueryStats queryStats;

    spdlog::logger *logger = opt.console.get();
    std::string dbg_file(opt.prefix + mantis::CQF_FILE);
    std::string sample_file(opt.prefix + mantis::SAMPLEID_FILE);

    std::vector<std::string> sampleNames = loadSampleFile(sample_file);
    queryStats.numSamples = sampleNames.size();
    logger->info("Number of experiments: {}", queryStats.numSamples);

    logger->info("Loading cqf...");
    CQF<KeyObject> cqf(dbg_file, CQF_FREAD);
    auto indexK = cqf.keybits() / 2;
    logger->info("Done loading cqf. k is {}", indexK);

    logger->info("Loading color classes...");
    MSTQuery mstQuery(opt.prefix, indexK, queryK, queryStats.numSamples, logger);
    logger->info("Done Loading color classes. Total # of color classes is {}",
                 mstQuery.parentbv.size() - 1);

    logger->info("Querying colored dbg.");
    std::ofstream opfile(opt.output);
    LRUCacheMap cache_lru(100000);
    RankScores rs(1);
    auto numOfQueries = mstQuery.parseKmers(opt.query_file, indexK);
    logger->info("Done parsing kmers in the query file.");
    mstQuery.findSamples(cqf, cache_lru, &rs, queryStats);
    logger->info("Done first round over the query file and finding samples for index k: {}", cqf.keybits()/2);
    if (opt.use_json) {
        output_results_json(opt.query_file, mstQuery, opfile, sampleNames, queryStats, numOfQueries);
    } else {
        output_results(opt.query_file, mstQuery, opfile, sampleNames, queryStats);
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