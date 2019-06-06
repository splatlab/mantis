//
// Created by Fatemeh Almodaresi on 2019-02-01.
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

struct Mc_stats {
    uint64_t nodeCnt = 1;
    uint64_t min_dist = -1; // infinity
    uint64_t color = -1; // infinity
};

struct workItem {
    dna::canonical_kmer node;
    colorIdType colorId;
    uint64_t kmerHash;
    uint64_t uniqueIdx;

    workItem(dna::canonical_kmer n, colorIdType c, uint64_t keyBits, uint64_t eqidx) : node(n), colorId(c),
    uniqueIdx(eqidx) {
        kmerHash = hash_64(n.val, BITMASK(keyBits));
//        std::cerr << "newworkItem:" << std::string(n) << " " << n.val << " " << kmerHash << "\n";
    }

    // Required to be able to use it as a key in set
    bool operator<(const workItem &item2) const {
        return kmerHash < item2.kmerHash;
    }
};

class Stat {
public:
    Stat(CQF<KeyObject>& cqfIn, uint64_t num_samples,
         spdlog::logger *logger): cqf(cqfIn), it(cqf.begin()) {
        k = cqf.keybits()/2;
        sdsl::util::assign(visited, sdsl::bit_vector(cqf.numslots(),0));
        std::cerr << "visited size: " << visited.size() << "\n";
    }
    Stat(CQF<KeyObject>& cqfIn, MSTQuery* mstQueryIn, uint64_t num_samples,
         spdlog::logger *logger): cqf(cqfIn), mstQuery(mstQueryIn), it(cqf.begin()) {
        cache_lru = new LRUCacheMap(100000);
        k = cqf.keybits()/2;
        oneCnt.resize((num_samples*(num_samples+1))/2);
        std::cout << "Total Eqs: " << mstQuery->parentbv.size() << "\n";
        sdsl::util::assign(visited, sdsl::bit_vector(cqf.numslots(),0));
    }
    std::vector<uint64_t> oneCnt;

    void operator++();
    Mc_stats operator*();
    bool done();

    std::vector<uint64_t> queryColor();
    void increaseCounter(uint64_t idx, uint64_t cnt);
    uint64_t getKey();
    uint64_t getColor();
private:
    uint32_t k;
    std::queue<workItem> work;
    std::unordered_set<uint64_t> visitedKeys;
    sdsl::bit_vector visited;
    MSTQuery* mstQuery;
    CQF<KeyObject>& cqf;
    CQF<KeyObject>::Iterator it;
    uint64_t num_samples;
    uint64_t kmerCntr = 0;
    LRUCacheMap* cache_lru;
    nonstd::optional<uint64_t> toDecode{nonstd::nullopt};
    QueryStats queryStats;
    std::set<workItem> neighbors(workItem n);
    bool exists(dna::canonical_kmer e, uint64_t &eqid, uint64_t &eqidx);
    uint64_t globalFuckingCounter = 0;
};

#endif //MANTIS_STAT_H
