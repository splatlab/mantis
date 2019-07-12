//
// Created by Fatemeh Almodaresi on 2019-07-12.
//

#ifndef MANTIS_SSKIQUERY_H
#define MANTIS_SSKIQUERY_H

#include <fstream>
#include <utility>
#include <spdlog/logger.h>

#include "sdsl/bit_vectors.hpp"
#include "canonicalKmerIterator.h"
#include "BooPHF.h"
#include "mantisconfig.hpp"

typedef boomphf::SingleHashFunctor<uint64_t> hasher_t;
typedef boomphf::mphf<uint64_t, hasher_t> boophf_t;


class SskiQuery {
public:
    bool queryKmer(CanonicalKmer kmer);
    SskiQuery(uint32_t kin, std::string prefix, spdlog::logger* c, bool tmpFlag)
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

    void parseKmers(std::string read, uint64_t kmer_size);

private:
    uint32_t k{23};
    sdsl::int_vector<2> contigSeq;
    sdsl::int_vector<3> prefixArr;
    std::unique_ptr<boophf_t> bphfPtr{nullptr};
    boophf_t* bphf;
    spdlog::logger* console;
    Sizes sizes;
    std::unordered_map<uint64_t, uint64_t> kmer2cidMap;

    uint64_t  binarySearch(CanonicalKmer kmer);
    uint64_t searchBucket(CanonicalKmer kmer, uint64_t idx);

};


#endif //MANTIS_SSKIQUERY_H
