//
// Created by Fatemeh Almodaresi on 7/20/18.
//

#include "MSF.h"
#include <unordered_set>
#include <random>

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

    MSFQuery(uint64_t numSamplesIn): numSamples(numSamplesIn) {
            numWrds = (uint64_t)std::ceil((double)numSamples/64.0);
    }

    void loadIdx(std::string indexDir) {
        sdsl::load_from_file(parentbv, indexDir+"/parents.bv");
        sdsl::load_from_file(deltabv, indexDir+"/deltas.bv");
        sdsl::load_from_file(bbv, indexDir+"/boundary.bv");
        sbbv = sdsl::bit_vector::select_1_type(&bbv);
        zero = parentbv.size()-1; // maximum color id which
        std::cerr << "Loaded the new color class index\n";
        std::cerr << "--> parent size: " << parentbv.size() << "\n"
                  << "--> delta size: " << deltabv.size() << "\n"
                  << "--> boundary size: " << bbv.size() << "\n";
        }

    std::vector<uint64_t> buildColor(uint64_t eqid, eqvec &bvs) {
        std::vector<uint32_t> flips(numSamples);
        uint64_t i{eqid}, from{0}, to{0};
        while (parentbv[i] != i) {
            //std::cerr << "c" << i << " p" << parentbv[i] << "\n";
            if (i > 0)
                from = sbbv(i)+1;
            else
                from = 0;
            to = sbbv(i+1);
            /*std::cerr << "delta list: \n";
            auto res = getDeltaList(bvs, parentbv[i], i, numSamples, numWrds);
            for (auto v : res) {
                std::cerr << v << " ";
            }
            std::cerr << "\nours from " << from  << " to " << to << ":\n";*/
            for (auto j = from; j <= to; j++) {
                //std::cerr << deltabv[j] << " ";
                flips[deltabv[j]] ^= 0x01;
            }
            //std::cerr <<"\n";
            i = parentbv[i];
        }
        if (i != zero) {
            //std::cerr << "root not zero\n";
            if (i > 0)
                from = sbbv(i)+1;
            to = sbbv(i+1);
            for (auto j = from; j <= to; j++) {
                //std::cerr << deltabv[j] << "\n";
                flips[deltabv[j]] ^= 0x01;
            }
        }
        /*else {
            std::cerr <<"root is zero\n";
        }*/
        std::vector<uint64_t> eq(numWrds);
        uint64_t one = 1;
        for (i = 0; i < numSamples; i++) {
            if (flips[i]) {
                uint64_t idx = i / 64;
                //std::cerr << "set " << i << " " << idx << " " << i % 64 << "\n";
                eq[idx] = eq[idx] | (one << (i % 64));
            }
        }
        return eq;
    }
};

uint32_t recursiveSteps(uint32_t idx, sdsl::int_vector<> &parentbv,
                        std::vector<uint32_t> &steps,
                        uint32_t zero) {
    if (idx == zero)
        return steps[idx]; // 0
    if (steps[idx] != 0)
        return steps[idx];
    if (parentbv[idx] == idx) {
        steps[idx]=1; // to retrieve the representative
        return steps[idx];
    }
    steps[idx] = recursiveSteps(parentbv[idx], parentbv, steps, zero)+1;
    return steps[idx];
}

struct Opts {
    std::string indexDir;
    uint64_t numSamples;
    std::string eqlistfile;
};

int main(int argc, char *argv[]) {

    using namespace clipp;
    enum class mode {
        validate, steps, query, help
    };
    mode selected = mode::help;
    Opts opt;

    auto validate_mode = (
            command("validate").set(selected, mode::validate),
                    required("-i", "--indexDir") &
                    value("index_dir", opt.indexDir) % "Directory containing index files.",
                    required("-eq", "--eqCls-lst-file") &
                    value("eqCls_list_filename", opt.eqlistfile) % "File containing list of equivalence (color) classes.",
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

    auto cli = (
            (validate_mode | steps_mode | command("help").set(selected, mode::help)
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

    std::string indexDir = opt.indexDir;
    std::string eqlistfile = opt.eqlistfile;
    uint64_t numSamples = opt.numSamples;
    uint64_t numWrds = (uint64_t)std::ceil((double)numSamples/64.0);
    MSFQuery query(numSamples);
    query.loadIdx(indexDir);
    uint64_t eqCount = query.parentbv.size()-1;
    std::cerr << "total # of equivalence classes is : " << eqCount << "\n";

    if (selected == mode::validate) {
        eqvec bvs;
        loadEqs(eqlistfile, bvs);
        // choose 1M random eq classes
        /*std::unordered_set<uint64_t> ids;
        while (ids.size() < 1000000) {
            ids.insert(rand() % eqCount);
        }*/

        uint64_t cntr{0};
        for (uint64_t idx = 0; idx < eqCount; idx++) {
            std::vector<uint64_t> newEq = query.buildColor(idx, bvs);
            std::vector<uint64_t> oldEq(numWrds);
            buildColor(bvs, oldEq, idx, numSamples);

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
    } else if (selected == mode::steps) {
        uint32_t zero = query.parentbv.size()-1;
        std::vector<uint32_t> steps(eqCount, 0);

        for (uint64_t idx = 0; idx < eqCount; idx++) {
            std::cout << query.parentbv[idx] << "\t"
                      << idx << "\t"
                      << recursiveSteps(idx, query.parentbv, steps, zero) << "\n";
        }
    }
}