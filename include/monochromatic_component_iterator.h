//
// Created by Fatemeh Almodaresi on 6/4/18.
//

#ifndef MANTIS_MONOCHROME_ITERATOR_H
#define MANTIS_MONOCHROME_ITERATOR_H

#include<set>
#include<queue>
#include <string>
#include <map>

#include <iostream>
#include <fstream>
#include <unordered_map>
#include <set>
#include <unordered_set>
#include <vector>

#include <inttypes.h>

#include "sparsepp/spp.h"
#include "tsl/sparse_map.h"
#include "sdsl/bit_vectors.hpp"
#include "bitvector.h"
#include "cqf.h"
#include "hashutil.h"
#include "common_types.h"

namespace dna {

/////////////// bases /////////////////
    enum base {
        C = 0, A = 1, T = 2, G = 3
    };

    base operator-(base b); // return the complementary base
    extern const base bases[4];
    extern const std::map<char, base> base_from_char;
    extern const std::map<base, char> base_to_char;

///////////// kmers /////////////////////
    class kmer {
    public:
        int len;
        uint64_t val;

        kmer(void);

        kmer(base b);

        kmer(int l, uint64_t v);

        kmer(std::string s);

        // Convert to string
        operator std::string() const;
    };

    bool operator<(kmer a, kmer b);

    bool operator==(kmer a, kmer b);

    bool operator!=(kmer a, kmer b);

// Return the reverse complement of k
    kmer operator-(kmer k);

    kmer canonicalize(kmer k);

// Return the kmer of length |a| that results from shifting b into a
// from the right
    kmer operator<<(kmer a, kmer b);

// Return the kmer of length |b| that results from shifting a into b
// from the left
    kmer operator>>(kmer a, kmer b);

// Append two kmers
    kmer operator+(kmer a, kmer b);

    kmer suffix(kmer k, int len);

    kmer prefix(kmer k, int len);

// The purpose of this class is to enable us to declare containers
// as holding canonical kmers, e.g. set<canonical_kmer>.  Then all
// inserts/queries/etc will automatically canonicalize their
// arguments.
    class canonical_kmer : public kmer {
    public:
        canonical_kmer(void);

        canonical_kmer(base b);

        canonical_kmer(int l, uint64_t v);

        canonical_kmer(std::string s);

        canonical_kmer(kmer k);
    };
}

struct Mc_stats {
    uint64_t nodeCnt = 0;
    uint64_t min_dist = -1; // infinity
};

typedef dna::canonical_kmer edge; // k-mer
typedef dna::canonical_kmer node; // (k-1)-mer

struct hash128 {
    uint64_t operator()(const __uint128_t &val128) const {
        __uint128_t val = val128;
        // Using the same seed as we use in k-mer hashing.
        return HashUtil::MurmurHash64A((void *) &val, sizeof(__uint128_t),
                                       2038074743);
    }
};

struct bvHash128 {
    __uint128_t operator()(const sdsl::bit_vector&bv) const {
        return HashUtil::MurmurHash128A((void *) bv.data(),
                                        bv.capacity()/8, 2038074743,
                                        2038074751);
    }
};

class monochromatic_component_iterator {
public:
    class work_item {
    public:
        node curr;
        uint64_t idx;
        uint64_t colorid;

        work_item(node currin, uint64_t idxin, uint64_t coloridin) : curr(currin), idx(idxin), colorid(coloridin) {}

        bool operator<(const work_item &item2) const {
            return (*this).curr < item2.curr;
        }
    };

    bool done();

    void operator++(void);

    Mc_stats operator*(void);

//monochromatic_component_iterator(const CQF<KeyObject> *g);
    monochromatic_component_iterator(const CQF<KeyObject> *g,
                                     BitVectorRRR &bvin,
                                     uint64_t num_samplesin = 2586);

    void neighborDist(uint64_t cntrr);
    void uniqNeighborDist(uint64_t num_samples);

    uint64_t cntr = 0;
    std::vector<uint64_t> withMax0;
    //spp::sparse_hash_map<__uint128_t, uint64_t, hash128> eqclass_map;
    spp::sparse_hash_map<sdsl::bit_vector, uint64_t, bvHash128> eqclass_map;


private:

    uint32_t k;
    std::queue<work_item> work;
    std::unordered_set<uint64_t> visitedKeys;
    const CQF<KeyObject> *cqf;
    CQF<KeyObject>::Iterator it;
    BitVectorRRR &bv;
    uint64_t num_samples;
    sdsl::bit_vector visited;

    bool exists(edge e, uint64_t &idx, uint64_t &eqid);

    std::set<monochromatic_component_iterator::work_item> neighbors(monochromatic_component_iterator::work_item n);

    work_item front(std::queue<work_item> &w);


    uint64_t manhattanDist(uint64_t eq1, uint64_t eq2);

    __uint128_t manhattanDistBvHash(uint64_t eq1, uint64_t eq2,
                                    uint64_t num_samples);
    void manhattanDistBvHash(uint64_t eq1, uint64_t eq2,
                                    sdsl::bit_vector& dist,
                                    uint64_t num_samples);


};


#endif //MANTIS_MONOCHROME_ITERATOR_H
