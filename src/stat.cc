//
// Created by Fatemeh Almodaresi on 2018-11-15.
//

#include <set>
#include "stat.h"
#include "ProgOpts.h"

void Stat::operator++(void) {

    if (it.done()) return; // don't cross the bound (undefined behaviour)
    ++it;
    kmerCntr++;
    //std::cerr << kmerCntr << "\r";
    if (kmerCntr % 1000000 == 0)
        std::cerr << "\rvisited " << kmerCntr << " kmers";
    while (!it.done()) {
        if (visitedKeys.find(it.get_cur_hash().key) != visitedKeys.end()) {
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
    workItem cur = {root, static_cast<colorIdType>((*it).count - 1), cqf.keybits()};
    work.push(cur);
    res.color = cur.colorId;
    visitedKeys.insert(cur.kmerHash);
    while (!work.empty()) {
        workItem w = work.front();
        work.pop();
        for (auto &neighbor : neighbors(w)) {
            //std::cerr << "nei";
            if (visitedKeys.find(neighbor.kmerHash) == visitedKeys.end()) {
                //std::cerr << "here:";
                if (neighbor.colorId == w.colorId) {
                    //std::cerr << "pup ";
                    work.push(neighbor);
                    res.nodeCnt++;
                    visitedKeys.insert(neighbor.kmerHash);
                }
            }
        }
    }
    return res;
}


std::set<workItem> Stat::neighbors(workItem n) {
    std::set<workItem> result;
    for (const auto b : dna::bases) {
        uint64_t eqid = 0;
        if (exists(n.node << b, eqid)) {
            result.insert(workItem(n.node << b, eqid, cqf.keybits()));
        }
        if (exists(b >> n.node, eqid)) {
            result.insert(workItem(b >> n.node, eqid, cqf.keybits()));
        }
    }
    return result;
}

bool Stat::exists(dna::canonical_kmer e, uint64_t &eqid) {
    KeyObject key(e.val, 0, 0);
    auto eqidtmp = cqf.query(key, 0);
    //std::cerr << "c" << eqidtmp;
    if (eqidtmp) {
        eqid = eqidtmp - 1;
        return true;
    }
    return false;
}

bool Stat::done() {return it.done();}

int stats_main(StatsOpts& sopt) {
    spdlog::logger *logger = sopt.console.get();
    std::string cqf_file = sopt.prefix + mantis::CQF_FILE;
    CQF<KeyObject> cqf(cqf_file, readmode::CQF_FREAD);
    Stat stats(cqf, sopt.numSamples, logger);
    if (sopt.type == "mono") {
        while (!stats.done()) {
            auto res = *stats;
            std::cout << res.color << "\t" << res.nodeCnt << "\n";
            ++stats;
        }
    }
    if (sopt.type == "cc_density") {
        uint32_t k = cqf.keybits()/2;
        MSTQuery mstQuery(sopt.prefix, k, k, sopt.numSamples, logger);
        Stat stats(cqf, &mstQuery, sopt.numSamples, logger);
        while (!stats.done()) {
            auto res = *stats;
            std::cout << res.color << "\t" << res.nodeCnt << "\n";
            ++stats;
        }
    }
}