#include <utility>

//
// Created by Fatemeh Almodaresi on 2019-03-07.
//

#ifndef MANTIS_MANTISSSKI_H
#define MANTIS_MANTISSSKI_H

#include "sdsl/bit_vectors.hpp"
#include "canonicalKmerIterator.h"
#include <climits>
#include "combine_kmer.h"
#include "string_view.hpp"
#include "BooPHF.h"


typedef boomphf::SingleHashFunctor<uint64_t> hasher_t;
typedef boomphf::mphf<uint64_t, hasher_t> boophf_t;

/**
 * Assumes the contig input is sorted lexicographically.
 * how I have sorted the output fasta file of bcalm in bash:
 * awk 'NR%2{printf "%s\t",$0;next;}1' unitig_fasta_file > tmp
 * sort --parallel=16 -t$'\t' -k2 tmp
 * sed 's/\t/\n/g' sorted_tmp > sorted_unitig_fasta_file
 */
class MantisSski {
// unitigVec
// startIdxVec
// prefixArr (3 bits)
// MPHf
public:
    MantisSski(uint32_t kin, std::string outd, spdlog::logger* c)
    : k(kin), outdir(std::move(outd)), console(c) {
        CanonicalKmer::k(kin);
        sizes = Sizes();
    }

    void buildUnitigVec(uint32_t numThreads, std::string rfile);
    void buildPrefixArr();
    void buildMPHF(uint32_t numThreads);

    bool queryKmer(CanonicalKmer kmer);

    void testBooPHF(uint32_t numThreads) {
        console->info("building boo ...");
        ContigKmerIterator kb(&contigSeq, sizes, static_cast<uint8_t>(k), 0);
        ContigKmerIterator ke(&contigSeq, sizes, static_cast<uint8_t>(k), contigSeq.size() - k + 1);
        std::vector<uint64_t> kmers;
        kmers.reserve(nkeys);
        while(kb != ke) {
            kmers.push_back(*kb);
            ++kb;
        }
        sort( kmers.begin(), kmers.end() );
        kmers.erase( unique( kmers.begin(), kmers.end() ), kmers.end() );
        console->info("done building kmers vector: {} = {}", nkeys, kmers.size());
        auto keyIt = boomphf::range(kmers.begin(), kmers.end());
        //std::cerr << " initialized the keyIt\n";
        bphf =
                new boophf_t(nkeys, keyIt, numThreads, 3.5); // keys.size(), keys, 16);
        console->info("mphf size = {}", (bphf->totalBitSize() / 8) / std::pow(2, 20));
    }



private:
    uint64_t contigAccumLength{0};
    uint64_t contigCnt{0};
    uint64_t nkeys{0};
    std::string outdir;
    uint32_t k{23};
    sdsl::int_vector<2> contigSeq;
    sdsl::int_vector<> contigStartIdx;
    sdsl::int_vector<3> prefixArr;
    boophf_t* bphf;
    spdlog::logger* console;
    Sizes sizes;
    void encodeSeq(sdsl::int_vector<2>& seqVec, size_t offset,
            stx::string_view str);

    bool isPalyndrome(std::string& st) {
        bool palyn = true;
        uint64_t i{0}, j{st.size()-1};
        while (palyn and i < j) {
            char r;
            switch (st[j]) {
                case 'A':
                    r = 'T';break;
                case 'C':
                    r = 'G';break;
                case 'G':
                    r = 'C';break;
                default:
                    r = 'A';
            }
            palyn = r == st[i];
            i++;
            j--;
        }
        return palyn;
    }
};



#endif //MANTIS_MANTISSSKI_H
