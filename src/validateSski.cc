//
// Created by Fatemeh Almodaresi on 2019-05-15.
//
#include "mantisconfig.hpp"
#include "sdsl/bit_vectors.hpp"
#include "gqf_cpp.h"
#include "canonicalKmer.h"
#include "canonicalKmerIterator.h"

int main(int argc, char* argv[]) {
    std::string prefix = argv[1];
    std::string kmer_file = argv[2];
    std::string cqf_file = prefix + "/" + mantis::CQF_FILE;
    std::string contigseq_file = prefix + "/" + mantis::SEQ_FILE;
    CQF<KeyObject> cqf(cqf_file, CQF_FREAD);
    uint32_t k = cqf.keybits() / 2;
    std::cerr << "cqf.numslots: " <<  cqf.numslots() << "\n";
    std::cerr << "cqf.dist_elts: " << cqf.dist_elts() << "\n";
    std::ofstream out(kmer_file, std::ios::out);
    sdsl::int_vector<2> contigSeq;
    sdsl::load_from_file(contigSeq, contigseq_file);
    Sizes sizes;
    uint64_t cntr{0};
    ContigKmerIterator kb(&contigSeq, sizes, static_cast<uint8_t>(k), 0);
    ContigKmerIterator ke(&contigSeq, sizes, static_cast<uint8_t>(k), contigSeq.size() - k + 1);
    while(kb != ke) {
        KeyObject key(*kb, 0, 0);
        auto eqidtmp = cqf.query(key, QF_NO_LOCK /*QF_KEY_IS_HASH | QF_NO_LOCK*/);
        if (!eqidtmp) {
            out << (*kb) << " f\n";
        } else {
            out << (*kb) << " nf\n";
        }
        ++kb;
        cntr++;
        if (cntr % 1000 == 0) std::cerr << "\r" << cntr/1000 << "k kmers processed";
    }

   /* CQF<KeyObject>::Iterator it(cqf.begin());
    while (!it.done()) {
        KeyObject keyobj = *it;
        dna::canonical_kmer kmer(static_cast<int>(k), keyobj.key);
        out << std::string(kmer) << "\n";
        ++it;
    }*/
    out.close();

}