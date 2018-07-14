/* -*- Mode: C++; c-basic-offset: 4; indent-tabs-mode: nil -*- */
#ifndef _HASHUTIL_H_
#define _HASHUTIL_H_

#include <sys/types.h>
#include <string>
#include <stdlib.h>
#include <stdint.h>

class HashUtil {
	public:

		// MurmurHash2
		static uint32_t MurmurHash(const void *buf, size_t length, uint32_t seed =
															 0);
		static uint32_t MurmurHash(const std::string &s, uint32_t seed = 0);
		static uint64_t MurmurHash64B ( const void * key, int len, unsigned int
																		seed );
		static uint64_t MurmurHash64A ( const void * key, int len, unsigned int
																		seed );
		static __uint128_t MurmurHash128A ( const void * key, int len, unsigned
																				int seed, unsigned int seed2 );

		// AES hash
		static uint64_t AES_HASH(uint64_t x);

		static uint64_t hash_64(uint64_t key, uint64_t mask);
		static uint64_t hash_64i(uint64_t key, uint64_t mask);

	private:
		HashUtil();
};

#endif  // #ifndef _HASHUTIL_H_


