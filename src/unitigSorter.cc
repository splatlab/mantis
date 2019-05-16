//
// Created by Fatemeh Almodaresi on 2019-05-14.
//
#include <vector>
#include <algorithm>
#include <fstream>

#include "combine_canonicalKmer.hpp"
#include "FastxParser.hpp"
#include "spdlog/spdlog.h"

#include "sdsl/bit_vectors.hpp"
#include "canonicalKmerIterator.h"

struct Unitig {
    std::string name;
    std::string seq;
    uint64_t kmerWrd;
    Unitig(std::string &n, std::string &s, uint64_t k):name{n}, seq{s}, kmerWrd{k} {}
};

int main(int argc, char* argv[]) {
    auto console = spdlog::stdout_color_mt("mantis_console");
    std::string cfile{argv[1]};
    std::string ofile{argv[2]};
    uint64_t k{std::stoul(argv[3])};
    uint64_t numThreads{2};
    if (argc > 4) numThreads = std::stoul(argv[4]);
    CanonicalKmer::k(k);
    // go over all the contig files
    console->info("Reading the contig file ...");
    std::vector<std::string> ref_files = {cfile};
    //FIXME how to use one parser to parse the file twice
    fastx_parser::FastxParser<fastx_parser::ReadSeq> parser(ref_files, numThreads, 1);
    parser.start();
    auto rg = parser.getReadGroup();
    std::vector<Unitig> unitgsFirstKmerVec;
    unitgsFirstKmerVec.reserve(1000000);
    uint64_t unitigCntr{0};
    Sizes sizes;
    // read the reference sequences and encode them into the refseq int_vector
    while (parser.refill(rg)) {
        for (auto &rp : rg) {
//            std::string seq = rp.seq;
            uint64_t curSize = rp.seq.size();
            uint64_t startIdx{0}, added2bits{0}, cntr{0}, quotient, remaining;
            do {
                remaining = curSize % (sizes.unitigSliceSize + 1);
                quotient = curSize / (sizes.unitigSliceSize + 1);
                curSize = curSize - sizes.unitigSliceSize + k-1;
                if (quotient) {
                    added2bits = sizes.unitigSliceSize;
                } else {
                    added2bits = remaining;
                }
                std::string seq = rp.seq.substr(startIdx, added2bits);
                startIdx += (added2bits-k+1);
                CanonicalKmer kmer;
                kmer.fromStr(seq.substr(0,k));
                std::string name = rp.name + std::to_string(cntr);
                unitgsFirstKmerVec.emplace_back(name, seq, kmer.fwWord());
                cntr++;
            } while (quotient);
        }
    }
    parser.stop();
    console->info("Sorting the unitigs...");
    std::sort(unitgsFirstKmerVec.begin(), unitgsFirstKmerVec.end(),
            [] (const auto& lhs, const auto& rhs){
        return lhs.kmerWrd < rhs.kmerWrd;
    });
    console->info("Writing the unitigs down on disk...");
    std::ofstream out(ofile);
    unitigCntr = 0;
    for (auto & p : unitgsFirstKmerVec) {
        out << ">" << p.name << " " << p.kmerWrd << "\n";
        out << p.seq << "\n";
    }
}