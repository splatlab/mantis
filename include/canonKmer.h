//
// Created by Fatemeh Almodaresi on 8/21/18.
//

#ifndef MANTIS_CANONKMER_H
#define MANTIS_CANONKMER_H

#include <map>
#include <iostream>
#include <string>
#include <assert.h>

#define BITMASK(nbits) ((nbits) == 64 ? 0xffffffffffffffff : (1ULL << (nbits)) \
												- 1ULL)

namespace duplicated_dna {

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

#endif //MANTIS_CANONKMER_H
