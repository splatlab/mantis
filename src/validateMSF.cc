//
// Created by Fatemeh Almodaresi on 7/20/18.
//

#include <unordered_set>
#include <random>
#include "MSF.h"
#include "cqf.h"
#include "common_types.h"
#include "CLI/CLI.hpp"
#include "CLI/Timer.hpp"
#include "kmer.h"

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
        zero = parentbv.size() - 1; // maximum color id which
        std::cerr << "Loaded the new color class index\n";
        std::cerr << "--> parent size: " << parentbv.size() << "\n"
                  << "--> delta size: " << deltabv.size() << "\n"
                  << "--> boundary size: " << bbv.size() << "\n";
    }

    std::vector<uint64_t> buildColor(uint64_t eqid, bool all=true) {
        std::vector<uint32_t> flips(numSamples);
        uint64_t i{eqid}, from{0}, to{0};
        while (parentbv[i] != i) {
            if (i > 0)
                from = sbbv(i) + 1;
            else
                from = 0;
            to = sbbv(i + 1);
            for (auto j = from; j <= to; j++) {
                flips[deltabv[j]] ^= 0x01;
            }
            i = parentbv[i];
        }
        if (i != zero) {
            if (i > 0)
                from = sbbv(i) + 1;
            to = sbbv(i + 1);
            for (auto j = from; j <= to; j++) {
                flips[deltabv[j]] ^= 0x01;
            }
        }
        if (!all) { // return the indices of set bits
            std::vector<uint64_t> eq;
            eq.reserve(numWrds);
            uint64_t one = 1;
            for (i = 0; i < numSamples; i++) {
                if (flips[i]) {
                    eq.push_back(i);
                }
            }
            return eq;
        }
        std::vector<uint64_t> eq(numWrds);
        uint64_t one = 1;
        for (i = 0; i < numSamples; i++) {
            if (flips[i]) {
                uint64_t idx = i / 64;
                eq[idx] = eq[idx] | (one << (i % 64));
            }
        }
        return eq;
    }

};

mantis::QueryResult findSamples(const mantis::QuerySet& kmers,
                               CQF<KeyObject>& dbg, MSFQuery& msfQuery) {
    std::unordered_map<uint64_t, uint64_t> query_eqclass_map;
    for (auto k : kmers) {
        KeyObject key(k, 0, 0);
        uint64_t eqclass = dbg.query(key);
        if (eqclass)
            query_eqclass_map[eqclass] += 1;
    }

    std::unordered_map<uint64_t, uint64_t> sample_map;
    for (auto it = query_eqclass_map.begin(); it != query_eqclass_map.end();
         ++it) {
        auto eqclass_id = it->first-1;
        auto count = it->second;
        auto setbits = msfQuery.buildColor(eqclass_id, false);
        for (auto sb : setbits) {
            sample_map[sb] += count;
        }
    }
    return sample_map;
}

uint32_t recursiveSteps(uint32_t idx, sdsl::int_vector<> &parentbv,
                        std::vector<uint32_t> &steps,
                        uint32_t zero) {
    if (idx == zero)
        return steps[idx]; // 0
    if (steps[idx] != 0)
        return steps[idx];
    if (parentbv[idx] == idx) {
        steps[idx] = 1; // to retrieve the representative
        return steps[idx];
    }
    steps[idx] = recursiveSteps(parentbv[idx], parentbv, steps, zero) + 1;
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

    if (selected == mode::validate) {
        eqvec bvs;
        loadEqs(opt.eqlistfile, bvs);
        uint64_t cntr{0};
        for (uint64_t idx = 0; idx < eqCount; idx++) {
            std::vector<uint64_t> newEq = msfQuery.buildColor(idx);
            std::vector<uint64_t> oldEq(numWrds);
            buildColor(bvs, oldEq, idx, opt.numSamples);
            if (newEq != oldEq) {
                std::cerr << "AAAAA! LOOOSER!!\n";
                std::cerr << cntr << ": index=" << idx << "\n";
                std::cerr << "new size: " << newEq.size() << " old size: " << oldEq.size() << "\n";
                for (auto k = 0; k < newEq.size(); k++) {
                    std::cerr << newEq[k] << " " << oldEq[k] << "\n";
                }
                std::exit(1);
            }
            cntr++;
            if (cntr % 10000000 == 0) {
                std::cerr << cntr << " eqs were the same\n";
            }
        }
        std::cerr << "WOOOOW! Validation passed\n";
    } else if (selected == mode::decodeAllEqs) {
        for (uint64_t idx = 0; idx < eqCount; idx++) {
            std::vector<uint64_t> newEq = msfQuery.buildColor(idx);
            if (idx % 10000000 == 0) {
                std::cerr << idx << " eqs decoded\n";
            }
        }
    } else if (selected == mode::steps) {
        uint32_t zero = msfQuery.parentbv.size() - 1;
        std::vector<uint32_t> steps(eqCount, 0);
        for (uint64_t idx = 0; idx < eqCount; idx++) {
            std::cout << msfQuery.parentbv[idx] << "\t"
                      << idx << "\t"
                      << recursiveSteps(idx, msfQuery.parentbv, steps, zero) << "\n";
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
        std::ofstream opfile(opt.outputfile);
        uint32_t cnt = 0;
        {
            CLI::AutoTimer timer{"Query time ", CLI::Timer::Big};
            for (auto &kmers : multi_kmers) {
                opfile << "seq" << cnt++ << '\t' << kmers.size() << '\n';
                mantis::QueryResult result = findSamples(kmers, dbg, msfQuery);
                for (auto it = result.begin(); it != result.end(); ++it) {
                    opfile << it->first/*cdbg.get_sample(it->first)*/ << '\t' << it->second << '\n';
                }
            }
            opfile.close();
        }
    }
}