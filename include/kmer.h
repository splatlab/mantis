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

#include "common_types.h"

#define BITMASK(nbits) ((nbits) == 64 ? 0xffffffffffffffff : (1ULL << (nbits)) \
												- 1ULL)
enum DNA_MAP {C, A, T, G};  // A=1, C=0, T=2, G=3

using namespace std;

class Kmer {
	public:
		static char map_int(uint8_t base);
		static uint8_t map_base(char base);
		static __int128_t str_to_int(std::string str);
		static std::string int_to_str(__int128_t kmer, uint64_t kmer_size);
		static int reverse_complement_base(int x);
		static __int128_t reverse_complement(__int128_t kmer, uint64_t kmer_size);
		static bool compare_kmers(__int128_t kmer, __int128_t kmer_rev);

		static mantis::QuerySets parse_kmers(const char *filename,
																				 uint64_t kmer_size, uint64_t&
																				 total_kmers);
		static std::string generate_random_string(uint64_t len);

	private:
		Kmer();
};
#endif
