//
// Created by Fatemeh Almodaresi on 2019-02-01.
//

#include <set>
#include <kmer.h>
#include "stat.h"
#include "ProgOpts.h"
#include "canonicalKmer.h"

void Stat::operator++(void) {

    if (it.done()) return; // don't cross the bound (undefined behaviour)
    ++it;
    kmerCntr++;
    //std::cerr << kmerCntr << "\r";
    if (kmerCntr % 1000000 == 0)
        std::cerr << "\rvisited " << kmerCntr << " kmers";
    while (!it.done()) {
        KeyObject key((*it).key, 0, 0);
        uint64_t eqidx = cqf.get_unique_index(key, 0);
        if (visited[eqidx]) {
            //if (visitedKeys.find(it.get_cur_hash().key) != visitedKeys.end()) {
            //if ((bool)(visited[it.iter.current])) {
            ++it;
            kmerCntr++;
            //std::cerr << kmerCntr << "\r";
            if (kmerCntr % 1000000 == 0)
                std::cerr << "\rvisited " << kmerCntr << " kmers";

        } else {
            break;
        }
    }
}


Mc_stats Stat::operator*(void) {
    if (!work.empty()) {
        std::cerr << "Throw Exception. The work queue should be empty at this point.\n";
        std::exit(1);
    }
    Mc_stats res;
    if (it.done()) return res;

    dna::canonical_kmer root(static_cast<int>(k), (*it).key);
    KeyObject key((*it).key, 0, 0);
    uint64_t eqidx = cqf.get_unique_index(key, 0);
    workItem cur = {root, static_cast<colorIdType>((*it).count - 1), cqf.keybits(), eqidx};
    work.push(cur);
    res.color = cur.colorId;
    while (!work.empty()) {
        workItem w = work.front();
        work.pop();
//        std::cerr << "v" << w.uniqueIdx << "\n";
        if (visited[w.uniqueIdx]) continue;
        visited[w.uniqueIdx] = 1;
        //if (visitedKeys.find(w.kmerHash) != visitedKeys.end()) continue;
        //visitedKeys.insert(w.kmerHash);

        for (auto &neighbor : neighbors(w)) {
            /* std::cerr << "n" << std::string(neighbor.node) << " " << neighbor.kmerHash
             << " " << neighbor.colorId << " "
             << neighbor.uniqueIdx << "\n";*/
            //std::cerr << "nei";
            if (visited[neighbor.uniqueIdx] == 0) {
                //if (visitedKeys.find(neighbor.kmerHash) == visitedKeys.end()) {
                /*std::cerr << globalFuckingCounter << "not visited\n";
                if (globalFuckingCounter > 50)
                    std::exit(1);
                globalFuckingCounter++;*/
                if (neighbor.colorId == w.colorId) {
                    //std::cerr << "pushed ";
                    work.push(neighbor);
                    res.nodeCnt++;
                    //visitedKeys.insert(neighbor.kmerHash);
                }
            }
            //std::cerr << "\n";
        }
    }
    return res;
}


std::set<workItem> Stat::neighbors(workItem n) {
    std::set<workItem> result;
    for (const auto b : dna::bases) {
        uint64_t eqid = 0, eqidx = 0;
        if (exists(n.node << b, eqid, eqidx)) {
            auto wi = workItem(n.node << b, eqid, cqf.keybits(), eqidx);
            //std::cerr<< "i" << wi.kmerHash << "\n";
            result.insert(wi);
        }
        if (exists(b >> n.node, eqid, eqidx)) {

            auto wi = workItem(b >> n.node, eqid, cqf.keybits(), eqidx);
            //std::cerr<< "i" << wi.kmerHash << "\n";
            result.insert(wi);
        }
    }
/*
    std::cerr << "nighbors:\n";
    for (auto res: result) {
        std::cerr << res.kmerHash << "\n";
    }
*/
    return result;
}

bool Stat::exists(dna::canonical_kmer e, uint64_t &eqid, uint64_t &eqidx) {
    KeyObject key(e.val, 0, 0);
    //std::cerr << "e:" << std::string(e) << " ";
    auto eqidtmp = cqf.query(key, 0);
    auto idx = cqf.get_unique_index(key, 0);

    //std::cerr << "c" << eqidtmp;
    if (eqidtmp) {
        eqid = eqidtmp - 1;
        eqidx = idx;
        //std::cerr << eqid << "\n";
        return true;
    }
    //std::cerr << "\n";
    return false;
}

bool Stat::done() { return it.done(); }

std::vector<uint64_t> Stat::queryColor() {
    colorIdType idx = static_cast<colorIdType>((*it).count - 1);
    std::vector<uint64_t> setbits;
    RankScores rs(1);
    nonstd::optional<uint64_t> dummy{nonstd::nullopt};

    if (cache_lru->contains(idx)) {
        setbits = (*cache_lru)[idx];//.get(eqclass_id);
        queryStats.cacheCntr++;
    } else {
        queryStats.noCacheCntr++;
        queryStats.trySample = (queryStats.noCacheCntr % 10 == 0);
        toDecode.reset();
        setbits = mstQuery->buildColor(idx, queryStats, cache_lru, &rs, toDecode);
        cache_lru->emplace(idx, setbits);
        if (queryStats.trySample and toDecode) {
            auto s = mstQuery->buildColor(*toDecode, queryStats, nullptr, nullptr, dummy);
            cache_lru->emplace(*toDecode, s);
        }
    }
    return setbits;
}

uint64_t Stat::getKey() { return (*it).key; }
uint64_t Stat::getColor() { return (*it).count - 1;}

void Stat::increaseCounter(uint64_t idx, uint64_t cnt) {
    uint64_t prefix = (idx * (idx + 1)) / 2 - 1;
    oneCnt[prefix + cnt]++;
}

int stats_main(StatsOpts &sopt) {
    spdlog::logger *logger = sopt.console.get();
    std::string cqf_file = sopt.prefix + mantis::CQF_FILE;
    CQF<KeyObject> cqf(cqf_file, readmode::CQF_FREAD);
    Stat stats(cqf, sopt.numSamples, logger);
    if (sopt.type == "mono") {
        std::unordered_map<uint64_t, std::vector<uint64_t>> mcc_freq;
        while (!stats.done()) {
            auto res = *stats;
            mcc_freq[res.nodeCnt].emplace_back(res.color);
            //std::cout << res.color << "\t" << res.nodeCnt << "\n";
            ++stats;
        }

        std::ofstream mcc_file(sopt.prefix + "/mcc_dist.out");
        mcc_file << "size_of_mcc\tcount\tcolor_dist\n";
        for (auto &kv : mcc_freq) {
            auto &mccs = kv.second;
            mcc_file << kv.first << "\t" << mccs.size() << "\t";
            std::sort(mccs.begin(), mccs.end(), [](uint64_t &v1, uint64_t &v2){
                return v1 < v2;
            });
            uint64_t prev{0}, color_cntr{0};
            for (auto c : mccs) {
                if (c == prev)
                    color_cntr++;
                else {
                    mcc_file << c << ":" << color_cntr << " ";
                    color_cntr = 0;
                }
                prev = c;
            }
        }
        mcc_file.close();
    }
    if (sopt.type == "color_dist") {
        std::unordered_map<uint64_t, uint64_t> colors;
        while(!stats.done()) {
            if (colors.find(stats.getColor()) == colors.end())
                colors[stats.getColor()] = 0;
            colors[stats.getColor()]+=1;
            ++stats;
        }

        std::ofstream dist_file(sopt.prefix + "/color_dist.out");
        dist_file << "color\tcount\n";
        for (auto &kv : colors) {
            dist_file << kv.first << "\t" << kv.second << "\n";
        }
        dist_file.close();
    }
    if (sopt.type == "cc_density") {
        uint32_t k = cqf.keybits() / 2;
        MSTQuery mstQuery(sopt.prefix, k, k, sopt.numSamples, logger);
        Stat stats(cqf, &mstQuery, sopt.numSamples, logger);
        while (!stats.done()) {
            std::vector<uint64_t> eq = stats.queryColor();
            uint64_t countOfOnes{1};
            //std::cerr << eq.size() << "\n";
            for (auto i : eq) {
                stats.increaseCounter(i, countOfOnes++);
            }
            ++stats;
        }
        for (auto cntr : stats.oneCnt) {
            std::cout << cntr << "\n";
        }
    }
    if (sopt.type == "jmerkmer") {
        uint32_t k = cqf.keybits() / 2;
        logger->info("k is {}", k);
        std::unordered_set<uint64_t> jmerMap;
        std::unordered_set<uint64_t> kmerMap;
        int j = (int) sopt.j, jmerCntr = j;
        logger->info("j is {}", j);
        uint8_t bases[] = {0, 1, 2, 3};
        uint64_t jmask = BITMASK(j);
        uint64_t jmaskExceptTheLowest2Bits = jmask - 3;
        uint64_t jmaskExceptTheHighest2Bits = BITMASK(j-2);
        uint64_t totalNumOfKmers = cqf.dist_elts(), kmerCntr=0;
        logger->info("Total number of kmers: {}", totalNumOfKmers);
        while (!stats.done()) {
            uint64_t originalKmer = stats.getKey();
            dna::canonical_kmer ckmer(k, originalKmer);
            uint64_t kmer = ckmer.val;
            if (kmerMap.find(kmer) != kmerMap.end()) continue;
            uint64_t jmer, dummy = BITMASK(64);
            jmerCntr = j;
//            std::cerr << "\nnew kmer\n";
            while (jmerCntr <= k) {
                bool branchIsFound = false;
                dna::canonical_kmer cjmer = dna::canonical_kmer(j, kmer);
                jmer = cjmer.val;
//                std::cerr << std::string(ckmer) << "\t" << std::string(cjmer) << "\n";
                if (jmerMap.find(jmer) == jmerMap.end()) {
                    for (uint64_t base : bases) { // replace 2 low order bits
                        uint64_t otherJmer = (jmer & jmaskExceptTheLowest2Bits) | base;
                        dna::canonical_kmer cother = dna::canonical_kmer(j, otherJmer);
//                        std::cerr << "other " << std::string(cother) << "\n";
                        if (jmerMap.find(dna::canonical_kmer(j, otherJmer).val) != jmerMap.end()) {
                            branchIsFound = true;
                            break;
                        }
                    }
                    if (!branchIsFound)
                        for (uint64_t base : bases) {
                            base <<= (2*(j-1)); // replace 2 high order bits
                            uint64_t otherJmer = (jmer & jmaskExceptTheHighest2Bits) | base;
                            if (jmerMap.find(dna::canonical_kmer(j, otherJmer).val) != jmerMap.end()) {
                                branchIsFound = true;
                                break;
                            }
                        }
                }
                if (branchIsFound) {
                    kmerMap.insert(ckmer.val);
                    break;
                } else {
                    jmerMap.insert(jmer);
                }
                kmer >>= 2;
                jmerCntr++;
            }
            ++stats;
            kmerCntr++;
            if (kmerCntr % 100000 == 0)
                std::cerr << "\r" << kmerCntr
                          << " confused kmers: " << kmerMap.size()
                          << " non-confusing jmers: " << jmerMap.size();
        }
        std::cerr << "\n";
        logger->info("total kmers observed: {}", kmerCntr);
        logger->info("total confused kmers: {}", kmerMap.size());
        logger->info("total non-confusing jmers: {}", jmerMap.size());
    }
}