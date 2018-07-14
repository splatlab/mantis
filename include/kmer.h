/*
 * =====================================================================================
 *
 *       Filename:  kmer.h
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  04/18/2016 05:06:51 PM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Prashant Pandey (), ppandey@cs.stonybrook.edu
 *   Organization:  Stony Brook University
 *
 * =====================================================================================
 */

#ifndef _KMER_H_
#define _KMER_H_

#include <stdio.h>
#include <string>

#include "hashutil.h"
#include "common_types.h"

constexpr const int K = 20;

#define BITMASK(nbits) ((nbits) == 64 ? 0xffffffffffffffff : (1ULL << (nbits)) \
												- 1ULL)
enum DNA_MAP {C, A, T, G};  // A=1, C=0, T=2, G=3

using namespace std;

class Kmer {
	public:
		static inline char map_int(uint8_t base);
		static inline uint8_t map_base(char base);
		static uint64_t str_to_int(string str);
		static string int_to_str(uint64_t kmer);
		static inline int reverse_complement_base(int x);
		static uint64_t reverse_complement(uint64_t kmer);
		static bool compare_kmers(uint64_t kmer, uint64_t kmer_rev);
		static inline unsigned __int128 word_reverse_complement(unsigned __int128 w);
		static inline int64_t word_reverse_complement(uint64_t w);
		static inline uint32_t word_reverse_complement(uint32_t w);
		static mantis::QuerySets parse_kmers(const char *filename,
																				 uint32_t seed, uint64_t range,
                                         uint64_t& total_kmers);
		static std::string generate_random_string(uint64_t len);

	private:
		Kmer();
};
#endif
