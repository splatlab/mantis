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

    std::vector<uint64_t> buildColor(uint64_t eqid) {
        std::vector<uint32_t> flips(numSamples);
        uint64_t i{eqid}, from{0}, to{0};
        while (parentbv[i] != i) {
            //std::cerr << i << " " << parentbv[i] << "\n";
            if (i > 0)
                from = sbbv(i)+1;
            to = sbbv(i+1);
            for (auto j = from; j <= to; j++) {
                //std::cerr << deltabv[j] << "\n";
                flips[deltabv[j]] ^= 0x01;
            }
            i = parentbv[i];
        }
        if (i != zero) {
            std::cerr << "root not zero\n";
            if (i > 0)
                from = sbbv(i)+1;
            to = sbbv(i+1);
            for (auto j = from; j <= to; j++) {
                //std::cerr << deltabv[j] << "\n";
                flips[deltabv[j]] ^= 0x01;
            }
        }
        else {
            std::cerr <<"root is zero\n";
        }
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


int main(int argc, char *argv[]) {
    std::string indexDir = argv[1];
    std::string eqlistfile = argv[2];
    uint64_t numSamples = std::stoull(argv[3]);
    uint64_t numWrds = (uint64_t)std::ceil((double)numSamples/64.0);
    MSFQuery query(numSamples);
    query.loadIdx(indexDir);
    uint64_t eqCount = query.parentbv.size()-1;
    std::cerr << "total # of equivalence classes is : " << eqCount << "\n";

    eqvec bvs;
    loadEqs(eqlistfile, bvs);
    // choose 1M random eq classes
    std::unordered_set<uint64_t> ids;
    while (ids.size() < 10) {
        ids.insert(rand() % eqCount);
    }

    uint64_t cntr{0};
    for (auto idx : ids) {
        std::vector<uint64_t> newEq = query.buildColor(idx);
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
    }
    std::cerr << "WOOOOW! Validation passed\n";
}