//
// Created by Fatemeh Almodaresi on 2018-10-13.
//

#include <assert.h>

#include "canonicalKmer.h"

#define BITMASK(nbits) ((nbits) == 64 ? 0xffffffffffffffff : (1ULL << (nbits)) - 1ULL)

namespace dna {

    /////////////// bases /////////////////
    base operator-(base b) {
        return (base) ((~((uint64_t) b)) & 0x3ULL);
    }

    const base bases[4] = {C, A, T, G};
    const std::map<char, base> base_from_char = {{'A', A},
                                                 {'C', C},
                                                 {'G', G},
                                                 {'T', T},
                                                 {'N', A}};
    const std::map<base, char> base_to_char = {{A, 'A'},
                                               {C, 'C'},
                                               {G, 'G'},
                                               {T, 'T'}};

    ///////////// kmers /////////////////////
    kmer::kmer(void) : len(0), val(0) {}

    kmer::kmer(base b) : len(1), val((uint64_t) b) {}

    kmer::kmer(int l, uint64_t v) : len(l), val(v & BITMASK(2 * l)) {
        assert(l <= 32);
    }

    static uint64_t string_to_kmer_val(std::string s) {
        uint64_t val = 0;
        for (auto c : s)
            val = (val << 2) | ((uint64_t) (base_from_char.at(c)));
        return val;
    }

    kmer::kmer(std::string s) : len(s.size()), val(string_to_kmer_val(s)) {
        assert(s.size() <= 32);
    }

    // Convert to string
    kmer::operator std::string() const {
        std::string s;
        for (auto i = 1; i < len + 1; i++)
            s = s + base_to_char.at((base) ((val >> (2 * (len - i))) & BITMASK(2)));
        return s;
    }

    bool operator<(kmer a, kmer b) {
        return a.len != b.len ? a.len < b.len : a.val < b.val;
    }

    bool operator==(kmer a, kmer b) {
        return a.len == b.len && a.val == b.val;
    }

    bool operator!=(kmer a, kmer b) {
        return !operator==(a, b);
    }

    // Return the reverse complement of k
    kmer operator-(kmer k) {
        uint64_t val = k.val;
        val =
                (val >> 32) |
                (val << 32);
        val =
                ((val >> 16) & 0x0000ffff0000ffff) |
                ((val << 16) & 0xffff0000ffff0000);
        val =
                ((val >> 8) & 0x00ff00ff00ff00ff) |
                ((val << 8) & 0xff00ff00ff00ff00);
        val =
                ((val >> 4) & 0x0f0f0f0f0f0f0f0f) |
                ((val << 4) & 0xf0f0f0f0f0f0f0f0);
        val =
                ((val >> 2) & 0x3333333333333333) |
                ((val << 2) & 0xcccccccccccccccc);
        val = ~val;
        val >>= 64 - 2 * k.len;
        return kmer(k.len, val);
    }

    // backwards from standard definition to match kmer.h definition
    kmer canonicalize(kmer k) {
        return -k < k ? k : -k;
    }

    // Return the kmer of length |a| that results from shifting b into a
    // from the right
    kmer operator<<(kmer a, kmer b) {
        uint64_t val = ((a.val << (2 * b.len)) | b.val) & BITMASK(2 * a.len);
        return kmer(a.len, val);
    }

    // Return the kmer of length |b| that results from shifting b into a
    // from the left
    kmer operator>>(kmer a, kmer b) {
        uint64_t val
                = ((b.val >> (2 * a.len)) | (a.val << (2 * (b.len - a.len))))
                  & BITMASK(2 * b.len);
        return kmer(b.len, val);
    }

    // Append two kmers
    kmer operator+(kmer a, kmer b) {
        int len = a.len + b.len;
        assert(len <= 32);
        uint64_t val = (a.val << (2 * b.len)) | b.val;
        return kmer(len, val);
    }

    kmer prefix(kmer k, int len) { return kmer(len, k.val >> (2 * (k.len - len))); }

    kmer suffix(kmer k, int len) { return kmer(len, k.val & BITMASK(2 * len)); }

    bool period_divides(kmer k, uint64_t periodicity) {
        static const uint64_t multipliers[33] =
                {
                        0,
                        0x5555555555555555, // 1
                        0x1111111111111111, // 2
                        0x1041041041041041, // 3
                        0x0101010101010101, // 4
                        0x1004010040100401, // 5
                        0x1001001001001001, // 6
                        0x0100040010004001, // 7
                        0x0001000100010001, // 8
                        0x0040001000040001, // 9
                        0x1000010000100001, // 10
                        0x0000100000400001, // 11
                        0x0001000001000001, // 12
                        0x0010000004000001, // 13
                        0x0100000010000001, // 14
                        0x1000000040000001, // 15
                        0x0000000100000001, // 16
                        0x0000000400000001, // 17
                        0x0000001000000001, // 18
                        0x0000004000000001, // 19
                        0x0000010000000001, // 20
                        0x0000040000000001, // 21
                        0x0000100000000001, // 22
                        0x0000400000000001, // 23
                        0x0001000000000001, // 24
                        0x0004000000000001, // 25
                        0x0010000000000001, // 26
                        0x0040000000000001, // 27
                        0x0100000000000001, // 28
                        0x0400000000000001, // 29
                        0x1000000000000001, // 30
                        0x4000000000000001, // 31
                        0x0000000000000001, // 32
                };
        uint64_t piece = k.val & BITMASK(2 * periodicity);
        piece = piece * multipliers[periodicity];
        piece = piece & BITMASK(2 * k.len);
        return piece == k.val;
    }

    uint64_t period(kmer k) {
        for (int i = 1; i <= k.len; i++) {
            if (period_divides(k, i))
                return i;
        }
        abort();
    }

    canonical_kmer::canonical_kmer(void) : kmer() {}

    canonical_kmer::canonical_kmer(base b) : kmer(canonicalize(kmer(b))) {}

    canonical_kmer::canonical_kmer(int l, uint64_t v)
            : kmer(canonicalize(kmer(l, v))) {}

    canonical_kmer::canonical_kmer(std::string s) : kmer(canonicalize(kmer(s))) {}

    canonical_kmer::canonical_kmer(kmer k) : kmer(canonicalize(k)) {}

}

