#include <utility>

//
// Created by Fatemeh Almodaresi on 2019-03-07.
//

#ifndef MANTIS_MANTISSSKI_H
#define MANTIS_MANTISSSKI_H

#include "sdsl/bit_vectors.hpp"
#include "canonicalKmerIterator.h"
#include <climits>
#include <fstream>

#include "combine_kmer.h"
#include "string_view.hpp"
#include "BooPHF.h"
#include "mantisconfig.hpp"


typedef boomphf::SingleHashFunctor<uint64_t> hasher_t;
typedef boomphf::mphf<uint64_t, hasher_t> boophf_t;

enum Cases {
    notFound,
    last2kmers,
    firstKmer
};

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

    MantisSski(uint32_t kin, std::string prefix, spdlog::logger* c, bool tmpFlag)
            : k(kin), console(c) {
        CanonicalKmer::k(kin);
        sizes = Sizes();
        console->info("Loading sequence vector...");
        sdsl::load_from_file(contigSeq, prefix + "/" + mantis::SEQ_FILE);

        {
//            CLI::AutoTimer timer{"Loading mphf table", CLI::Timer::Big};
            console->info("Loading mphf table...");
            std::string hfile = prefix + "/" + mantis::MPHF_FILE;
            std::ifstream hstream(hfile);
            bphfPtr.reset(new boophf_t);
            bphfPtr->load(hstream);
            hstream.close();
            bphf = bphfPtr.get();
        }
        console->info("Loading prefix array...");
        sdsl::load_from_file(prefixArr, prefix + "/" + mantis::PREFIX_FILE);
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

    uint64_t  binarySearch(CanonicalKmer kmer);
    bool searchBucket(CanonicalKmer kmer, uint64_t idx);
    uint64_t contigCnt{0};
    uint64_t nkeys{0};
    std::string outdir;
    uint32_t k{23};
    sdsl::int_vector<2> contigSeq;
    sdsl::int_vector<3> prefixArr;
    std::unique_ptr<boophf_t> bphfPtr{nullptr};
    boophf_t* bphf;
    spdlog::logger* console;
    Sizes sizes;
};



#endif //MANTIS_MANTISSSKI_H
