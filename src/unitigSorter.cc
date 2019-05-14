//
// Created by Fatemeh Almodaresi on 2019-05-14.
//
#include <vector>
#include <algorithm>
#include <fstream>

#include "combine_canonicalKmer.hpp"
#include "FastxParser.hpp"

struct Unitig {
    std::string name;
    std::string seq;
    uint64_t kmerWrd;
    Unitig(std::string &n, std::string &s, uint64_t k):name{n}, seq{s}, kmerWrd{k} {}
};

int main(int argc, char* argv[]) {

    std::string cfile{argv[1]};
    std::string ofile{argv[2]};
    uint64_t k{std::stoul(argv[3])};
    uint64_t numThreads{2};
    if (argc > 4) numThreads = std::stoul(argv[4]);
    CanonicalKmer::k(k);
    // go over all the contig files
    std::cerr << "Reading the contig file ...\n";
    std::vector<std::string> ref_files = {cfile};
    //FIXME how to use one parser to parse the file twice
    fastx_parser::FastxParser<fastx_parser::ReadSeq> parser(ref_files, numThreads, 1);
    parser.start();
    auto rg = parser.getReadGroup();
    std::vector<Unitig> unitgsFirstKmerVec;
    unitgsFirstKmerVec.reserve(1000000);
    uint64_t unitigCntr{0};
    // read the reference sequences and encode them into the refseq int_vector
    while (parser.refill(rg)) {
        for (auto &rp : rg) {
            std::string seq = rp.seq;
            CanonicalKmer kmer;
            kmer.fromStr(seq.substr(0,k));
            unitgsFirstKmerVec.emplace_back(rp.name, rp.seq, kmer.fwWord());
            unitigCntr++;
        }
    }
    parser.stop();
    std::cerr << "Sorting the unitigs..\n";
    std::sort(unitgsFirstKmerVec.begin(), unitgsFirstKmerVec.end(),
            [] (const auto& lhs, const auto& rhs){
        return lhs.kmerWrd < rhs.kmerWrd;
    });
    std::cerr << "Writing the unitigs down on disk..\n";
    std::ofstream out(ofile);
    unitigCntr = 0;
    for (auto & p : unitgsFirstKmerVec) {
        out << ">" << p.name << " " << p.kmerWrd << "\n";
        out << p.seq << "\n";
    }
}