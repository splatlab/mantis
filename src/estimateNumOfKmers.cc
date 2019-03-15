//
// Created by Fatemeh Almodaresi on 2018-10-15.
//

#include <iostream>
#include <unordered_set>
#include <string>
#include <fstream>

#include "canonicalKmer.h"
#include "gqf/hashutil.h"

#define BITMASK(nbits)  ((nbits) == 64 ? 0xffffffffffffffff : (1ULL << (nbits)) - 1ULL)

int main(int argc, char* argv[]) {
    if (argc < 2) {
        std::cerr << "Require two arguments, "
                     "At least list of ntCard files and also value of k (optional, default 31).\n";
        std::exit(1);
    }
    std::string inputFile = argv[1];
    uint64_t k = 31;
    if (argc == 3)
        k = std::stoull(argv[2]);
    std::ifstream inf(inputFile);
    std::unordered_set<uint64_t> kmerSet;
    std::string jellyfishFile;
    uint64_t cntr{0}, fileCntr{0};
    while (inf >> jellyfishFile) {
        std::ifstream jellyfish(jellyfishFile);
        std::string kmer;
        uint64_t cnt;
        while (jellyfish >> kmer >> cnt) {
            cntr++;
            dna::canonical_kmer ckmer(kmer);
            kmerSet.insert(hash_64(ckmer.val, BITMASK(2*k)));
        }
        jellyfish.close();
        std::cerr << "\r" << ++fileCntr << " files, " << cntr/1000000 << "M kmers, " << kmerSet.size() << " unique kmers";
    }
    inf.close();
    std::cerr << "\r\nObserved " << fileCntr << " files, " <<
    cntr/1000000 << "M kmers and " << kmerSet.size() << " unique kmers\n";
}