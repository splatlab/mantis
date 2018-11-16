//
// Created by Fatemeh Almodaresi on 2018-11-15.
//

#ifndef MANTIS_STAT_H
#define MANTIS_STAT_H

#include <vector>
#include <unordered_set>
#include <queue>

#include "gqf_cpp.h"
#include "lru/lru.hpp"
#include "canonicalKmer.h"
#include "mstQuery.h"
#include "gqf/hashutil.h"

typedef uint32_t colorIdType;

#define MAX_VALUE(nbits) ((1ULL << (nbits)) - 1)
#define BITMASK(nbits)                                    \
  ((nbits) == 64 ? 0xffffffffffffffff : MAX_VALUE(nbits))

struct Mc_stats {
    uint64_t nodeCnt = 1;
    uint64_t min_dist = -1; // infinity
    uint64_t color = -1; // infinity
};

struct workItem {
    dna::canonical_kmer node;
    colorIdType colorId;
    uint64_t kmerHash;

    workItem(dna::canonical_kmer n, colorIdType c, uint64_t keyBits) : node(n), colorId(c) {
        kmerHash = hash_64(n.val, BITMASK(keyBits));
    }

    // Required to be able to use it as a key in set
    bool operator<(const workItem &item2) const {
        return (*this).node < item2.node;
    }
};

class Stat {
public:
    Stat(CQF<KeyObject>& cqfIn, uint64_t num_samples,
         spdlog::logger *logger): cqf(cqfIn), it(cqf.begin()) {
        k = cqf.keybits()/2;
    }
    Stat(CQF<KeyObject>& cqfIn, MSTQuery* mstQueryIn, uint64_t num_samples,
         spdlog::logger *logger): cqf(cqfIn), mstQuery(mstQueryIn), it(cqf.begin()) {
        k = cqf.keybits()/2;
    }
    void operator++();
    Mc_stats operator*();
    bool done();

private:
    uint32_t k;
    std::queue<workItem> work;
    std::unordered_set<uint64_t> visitedKeys;
    MSTQuery* mstQuery;
    CQF<KeyObject>& cqf;
    CQF<KeyObject>::Iterator it;
    uint64_t num_samples;
    uint64_t kmerCntr = 0;

    std::set<workItem> neighbors(workItem n);
    bool exists(dna::canonical_kmer e, uint64_t &eqid);

};


#endif //MANTIS_STAT_H
