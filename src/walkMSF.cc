//
// Created by Fatemeh Almodaresi on 7/20/18.
//

#include <unordered_set>
#include <random>
#include <chrono>
#include "MSF.h"
#include "cqf.h"
#include "common_types.h"
#include "CLI/CLI.hpp"
#include "CLI/Timer.hpp"
#include "kmer.h"
#include "lrucache.hpp"
#include "hashutil.h"
#include "lru/lru.hpp"
#include "tsl/hopscotch_map.h"
#include "nonstd/optional.hpp"
#include <bitset>

struct QueryStats {
    uint32_t cnt = 0, cacheCntr = 0, noCacheCntr{0};
    uint64_t totSel{0};
    std::chrono::duration<double> selectTime{0};
    std::chrono::duration<double> flipTime{0};
    uint64_t totEqcls{0};
    uint64_t rootedNonZero{0};
    uint64_t nextCacheUpdate{10000};
    uint64_t globalQueryNum{0};
    std::vector<uint64_t> buffer;
    uint64_t numSamples{0};
  tsl::hopscotch_map<uint32_t, uint64_t> numOcc;
  bool trySample{false};
};


namespace mantis{
  namespace util {
    class int_hasher {
    public:
      size_t operator()(uint64_t i) const {
        return HashUtil::MurmurHash64A(reinterpret_cast<const void*>(&i), sizeof(i), 8675309);
      }
    };
  }
}

using LRUCacheMap =  LRU::Cache<uint64_t, std::vector<uint64_t>>;// cache::lru_cache <uint64_t, std::vector<uint64_t>, mantis::util::int_hasher>;

class RankScores {
public:
  RankScores(uint64_t nranks) {rs_.resize(nranks);}

  std::unordered_map<uint32_t, uint32_t>& operator[](uint32_t r) {
    if (r > maxRank_) {
      maxRank_ = std::min(r, static_cast<uint32_t>(rs_.size()-1));
    }
    return (r < rs_.size()) ? rs_[r] : rs_.back();
  }

  void clear() {
    for (auto& m : rs_){ m.clear(); }
    maxRank_ = 0;
  }

  uint32_t maxRank() const { return maxRank_; }

private:
  std::vector<std::unordered_map<uint32_t, uint32_t>> rs_;
  uint32_t maxRank_{0};
};


class MSFQuery {

private:
    uint64_t numSamples;
    uint64_t numWrds;
    uint32_t zero;
    sdsl::bit_vector bbv;

public:
    sdsl::int_vector<> parentbv;
    sdsl::int_vector<> deltabv;
    sdsl::bit_vector::select_1_type sbbv;

    MSFQuery(uint64_t numSamplesIn) : numSamples(numSamplesIn) {
        numWrds = (uint64_t) std::ceil((double) numSamples / 64.0);
    }

    void loadIdx(std::string indexDir) {
        sdsl::load_from_file(parentbv, indexDir + "/parents.bv");
        sdsl::load_from_file(deltabv, indexDir + "/deltas.bv");
        sdsl::load_from_file(bbv, indexDir + "/boundary.bv");
        sbbv = sdsl::bit_vector::select_1_type(&bbv);
        zero = 0; //parentbv.size() - 1; // maximum color id which
        std::cerr << "Loaded the new color class index\n";
        std::cerr << "--> parent size: " << parentbv.size() << "\n"
                  << "--> delta size: " << deltabv.size() << "\n"
                  << "--> boundary size: " << bbv.size() << "\n";
    }

    std::vector<uint64_t> buildColor(uint64_t eqid, QueryStats &queryStats,
                                     LRUCacheMap *lru_cache,
                                     RankScores* rs,
                                     nonstd::optional<uint64_t>& toDecode, // output param.  Also decode these
                                     bool all = true) {
        eqid++;
        (void)rs;
        std::vector<uint32_t> flips(numSamples);
        std::vector<uint32_t> xorflips(numSamples, 0);
        uint64_t i{eqid}, from{0}, to{0};
        int64_t height{0};
        auto& froms = queryStats.buffer;
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
              auto& occ = queryStats.numOcc[iparent];
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
            if (flips[i] ^  xorflips[i]) {
                uint64_t idx = i / 64;
                eq[idx] = eq[idx] | (one << (i % 64));
            }
        }
        return eq;
    }

};

mantis::QueryResult findSamples(const mantis::QuerySet &kmers,
                                CQF<KeyObject> &dbg, MSFQuery &msfQuery,
                                LRUCacheMap& lru_cache,
                                RankScores* rs,
                                QueryStats &queryStats) {
    std::unordered_map<uint64_t, uint64_t> query_eqclass_map;
    for (auto k : kmers) {
        KeyObject key(k, 0, 0);
        uint64_t eqclass = dbg.query(key);
        if (eqclass)
            query_eqclass_map[eqclass] += 1;
    }

    mantis::QueryResult sample_map(queryStats.numSamples,0);
    size_t numPerLevel = 10;
    nonstd::optional<uint64_t> toDecode{nonstd::nullopt};
    nonstd::optional<uint64_t> dummy{nonstd::nullopt};

    for (auto it = query_eqclass_map.begin(); it != query_eqclass_map.end(); ++it) {
        auto eqclass_id = it->first - 1;
        auto count = it->second;

        std::vector<uint64_t> setbits;
        if (lru_cache.contains(eqclass_id)) {
            setbits = lru_cache[eqclass_id];//.get(eqclass_id);
            queryStats.cacheCntr++;
        } else {
            queryStats.noCacheCntr++;
            toDecode.reset();
            dummy.reset();
            queryStats.trySample = (queryStats.noCacheCntr % 10 == 0);
            setbits = msfQuery.buildColor(eqclass_id, queryStats, &lru_cache, rs, toDecode, false);
            lru_cache.emplace(eqclass_id, setbits);
            if ((queryStats.trySample) and toDecode) {
              auto s = msfQuery.buildColor(*toDecode, queryStats, nullptr, nullptr, dummy, false);
              lru_cache.emplace(*toDecode, s);
            }
        }
        for (auto sb : setbits) {
            sample_map[sb] += count;
        }

        ++queryStats.globalQueryNum;
        /*
        if (queryStats.globalQueryNum > queryStats.nextCacheUpdate) {
          for (int64_t i = rs.maxRank(); i > 50; i-=50) {
            auto& m = rs[i];
            if (m.size() > numPerLevel) {
              std::vector<std::pair<uint32_t, uint32_t>> pairs;
              pairs.reserve(m.size());
              std::copy(m.begin(), m.end(), std::back_inserter(pairs));
              std::nth_element(pairs.begin(), pairs.begin()+numPerLevel, pairs.end(),
                               [](const std::pair<uint32_t, uint32_t> &a, const std::pair<uint32_t, uint32_t> &b) {
                                 return a.second > b.second;
                               });
              for (auto pit = pairs.begin(); pit != pairs.begin()+numPerLevel; ++pit) {
                if(!lru_cache.exists(pit->first)) {
                  auto v = msfQuery.buildColor(pit->first, queryStats, nullptr, nullptr, false);
                  lru_cache.put(pit->first, v);
                }
              }
            } else if(m.size() > 0){
              for (auto pit = m.begin(); pit != m.end(); ++pit) {
                if(!lru_cache.exists(pit->first)) {
                  auto v = msfQuery.buildColor(pit->first, queryStats, nullptr, nullptr, false);
                  lru_cache.put(pit->first, v);
                }
              }
            }

          }
          queryStats.nextCacheUpdate += 10000;
          rs.clear();
        }*/
    }

    return sample_map;
}

std::pair<uint32_t, uint32_t> recursiveSteps(uint32_t idx, sdsl::int_vector<> &parentbv,
                                             std::vector<std::pair<uint32_t, uint32_t>> &steps,
                                             uint32_t zero) {
    if (idx == zero) {
        steps[idx].second = zero;
        return steps[idx]; // 0
    }
    if (steps[idx].first != 0)
        return steps[idx];
    if (parentbv[idx] == idx) {
        steps[idx].first = 1; // to retrieve the representative
        steps[idx].second = idx;
        return steps[idx];
    }
    auto ret = recursiveSteps(parentbv[idx], parentbv, steps, zero);
    steps[idx].first = ret.first + 1;
    steps[idx].second = ret.second;
    return steps[idx];
}

struct Opts {
    std::string indexDir;
    uint64_t numSamples;
    std::string eqlistfile;
    std::string cqffile;
    std::string outputfile;
    std::string queryfile;
};

int main(int argc, char *argv[]) {
    ios_base::sync_with_stdio(false);
    using namespace clipp;
    enum class mode {
        validate, steps, decodeAllEqs, query, help
    };
    mode selected = mode::help;
    Opts opt;

    auto validate_mode = (
            command("validate").set(selected, mode::validate),
                    required("-i", "--indexDir") &
                    value("index_dir", opt.indexDir) % "Directory containing index files.",
                    required("-eq", "--eqCls-lst-file") &
                    value("eqCls_list_filename", opt.eqlistfile) %
                    "File containing list of equivalence (color) classes.",
                    required("-s", "--numSamples") &
                    value("numSamples", opt.numSamples) % "Total number of experiments (samples)."
    );

    auto steps_mode = (
            command("steps").set(selected, mode::steps),
                    required("-i", "--indexDir") &
                    value("index_dir", opt.indexDir) % "Directory containing index files.",
                    required("-s", "--numSamples") &
                    value("numSamples", opt.numSamples) % "Total number of experiments (samples)."
    );

    auto decodeAllEqs_mode = (
            command("decodeAllEqs").set(selected, mode::decodeAllEqs),
                    required("-i", "--indexDir") &
                    value("index_dir", opt.indexDir) % "Directory containing index files.",
                    required("-s", "--numSamples") &
                    value("numSamples", opt.numSamples) % "Total number of experiments (samples)."
    );

    auto query_mode = (
            command("query").set(selected, mode::query),
                    required("-i", "--indexDir") &
                    value("index_dir", opt.indexDir) % "Directory containing index files.",
                    required("-g", "--cqf") &
                    value("kmer-graph", opt.cqffile) % "cqf file containing the kmer mapping to color class ids.",
                    required("-q", "--queryFile") &
                    value("query-file", opt.queryfile) % "query file containing list of sequences.",
                    required("-o", "--outputFile") &
                    value("query-output-file", opt.outputfile) % "file to write query results.",
                    required("-s", "--numSamples") &
                    value("numSamples", opt.numSamples) % "Total number of experiments (samples)."
    );

    auto cli = (
            (validate_mode | steps_mode | decodeAllEqs_mode | query_mode | command("help").set(selected, mode::help)
            )
    );

    decltype(parse(argc, argv, cli)) res;
    try {
        res = parse(argc, argv, cli);
    } catch (std::exception &e) {
        std::cout << "\n\nParsing command line failed with exception: " << e.what() << "\n";
        std::cout << "\n\n";
        std::cout << make_man_page(cli, "MSF");
        return 1;
    }
    if (!res) {
        std::cerr << "Cannot parse the input arguments\n";
        std::exit(1);
    }
    if (selected == mode::help) {
        std::cerr << make_man_page(cli, "MSF");
        std::exit(1);
    }

    uint64_t numWrds = (uint64_t) std::ceil((double) opt.numSamples / 64.0);
    MSFQuery msfQuery(opt.numSamples);
    msfQuery.loadIdx(opt.indexDir);
    uint64_t eqCount = msfQuery.parentbv.size() - 1;
    std::cerr << "total # of equivalence classes is : " << eqCount << "\n";

    LRUCacheMap cache_lru(100000);
    QueryStats queryStats;
    queryStats.numSamples = opt.numSamples;

    if (selected == mode::validate) {
        eqvec bvs;
        loadEqs(opt.eqlistfile, bvs);
        uint64_t cntr{0};
        for (uint64_t idx = 0; idx < eqCount; idx++) {
          nonstd::optional<uint64_t> dummy{nonstd::nullopt};
            std::vector<uint64_t> newEq = msfQuery.buildColor(idx, queryStats, &cache_lru, nullptr, dummy);
            std::vector<uint64_t> oldEq(numWrds);
            buildColor(bvs, oldEq, idx, opt.numSamples);
            if (newEq != oldEq) {
                std::cerr << "AAAAA! LOOOSER!!\n";
                std::cerr << cntr << ": index=" << idx << "\n";
                std::cerr << "n ";
                for (auto k = 0; k < newEq.size(); k++) {
                    std::cerr << std::bitset<64>(newEq[k]);
                }
                std::cerr << "\no ";
                for (auto k = 0; k < oldEq.size(); k++) {
                    std::cerr << std::bitset<64>(oldEq[k]);
                }
                std::cerr << "\n";
                std::exit(1);
            }
            cntr++;
            if (cntr % 10000000 == 0) {
                std::cerr << cntr << " eqs were the same\n";
            }
        }
        std::cerr << "WOOOOW! Validation passed\n";
    } else if (selected == mode::decodeAllEqs) {
        std::random_device r;

        // Choose a random mean between 1 and 6
        std::default_random_engine e1(r());
        std::uniform_int_distribution<int> uniform_dist(0, eqCount - 1);
        for (uint64_t idx = 0; idx < 182169; idx++) {
          nonstd::optional<uint64_t> dummy{nonstd::nullopt};
            std::vector<uint64_t> newEq = msfQuery.buildColor(uniform_dist(e1),
                                                              queryStats,
                                                              &cache_lru,
                                                              nullptr,
                                                              dummy,
                                                              false);
            /*if (idx % 10000000 == 0) {
                std::cerr << idx << " eqs decoded\n";
            }*/
        }
        std::cerr << "cache was used " << queryStats.cacheCntr << " times " << queryStats.noCacheCntr << "\n";
        //std::cerr << "select time was " << queryStats.selectTime.count() << "s, flip time was "
        //         << queryStats.flipTime.count() << '\n';
        //std::cerr << "total selects = " << queryStats.totSel << ", time per select = "
        //          << queryStats.selectTime.count() / queryStats.totSel << '\n';
        std::cerr << "total # of queries = " << queryStats.totEqcls
                  << ", total # of queries rooted at a non-zero node = " << queryStats.rootedNonZero << "\n";

    } else if (selected == mode::steps) {
        uint32_t zero = msfQuery.parentbv.size() - 1;
        std::vector<std::pair<uint32_t, uint32_t>> steps(eqCount, {0, 0});
        for (uint64_t idx = 0; idx < eqCount; idx++) {
            auto row = recursiveSteps(idx, msfQuery.parentbv, steps, zero);
            std::cout << msfQuery.parentbv[idx] << "\t"
                      << idx << "\t"
                      << row.first << "\t" << row.second << "\n";
        }
    } else if (selected == mode::query) {
        CQF<KeyObject> dbg(opt.cqffile, false);
        std::cerr << "Done loading cqf.\n";
        uint32_t seed = 2038074743;
        uint64_t total_kmers = 0;
        // loading kmers
        mantis::QuerySets multi_kmers = Kmer::parse_kmers(opt.queryfile.c_str(),
                                                          seed,
                                                          dbg.range(),
                                                          20,
                                                          total_kmers);
        std::cerr << "Done loading query file : # of seqs: " << multi_kmers.size() << "\n";
        RankScores rs(1);
        std::ofstream opfile(opt.outputfile);
        {
            CLI::AutoTimer timer{"Query time ", CLI::Timer::Big};
            for (auto &kmers : multi_kmers) {
                opfile << "seq" << queryStats.cnt++ << '\t' << kmers.size() << '\n';
                mantis::QueryResult result = findSamples(kmers, dbg, msfQuery, cache_lru,
                                                         nullptr, queryStats);
                for (auto it = result.begin(); it != result.end(); ++it) {
                  if (*it > 0) {
                  auto i = std::distance(result.begin(), it);
                  opfile << i << '\t' << *it << '\n';
                  }
                  //opfile << it->first/*cdbg.get_sample(it->first)*/ << '\t' << it->second << '\n';
                }
            }
            opfile.close();
        }
        std::cerr << "cache was used " << queryStats.cacheCntr << " times " << queryStats.noCacheCntr << "\n";
        std::cerr << "select time was " << queryStats.selectTime.count() << "s, flip time was "
                  << queryStats.flipTime.count() << '\n';
        std::cerr << "total selects = " << queryStats.totSel << ", time per select = "
                  << queryStats.selectTime.count() / queryStats.totSel << '\n';
        std::cerr << "total # of queries = " << queryStats.totEqcls
                  << ", total # of queries rooted at a non-zero node = " << queryStats.rootedNonZero << "\n";
        /*for (auto &kv : queryStats.numOcc) {
            std::cout << kv.first << '\t' << kv.second << '\n';
        }*/
    }
}
