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
    if (argc < 3) {
        console->error("Require at least 3 arguments: input contig file, "
                      "output contig file, and k are required and numThreads is optional (Default=4).");
        std::exit(3);
    }
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
            uint64_t startIdx{0}, added2bits{0}, cntr{0}, quotient, remaining;
            std::string seq = rp.seq;
            // this is the ridiculous condition for taking care of a corner case in bcalm
            if (seq.size() > k and seq.substr(0, k) == seq.substr(seq.size() - k, k)) {
                seq = seq.substr(0, seq.size() - 1);
            }
            uint64_t curSize = seq.size();
//            if (rp.seq.size() > k and isPalyndrome(rp.seq))
//                std::cout << rp.seq << "\n";


            do {
                remaining = curSize % (sizes.unitigSliceSize + 1);
                quotient = curSize / (sizes.unitigSliceSize + 1);
                curSize = curSize - sizes.unitigSliceSize + k-1;
                if (quotient) {
                    added2bits = sizes.unitigSliceSize;
                } else {
                    added2bits = remaining;
                }
                std::string currseq = seq.substr(startIdx, added2bits);
                startIdx += (added2bits-k+1);
                CanonicalKmer kmer;
                kmer.fromStr(currseq.substr(0,k));
                std::string name = rp.name + std::to_string(cntr);
                unitgsFirstKmerVec.emplace_back(name, currseq, kmer.fwWord());
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