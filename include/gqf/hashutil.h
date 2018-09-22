/*
 * ============================================================================
 *
 *        Authors:  Prashant Pandey <ppandey@cs.stonybrook.edu>
 *                  Rob Johnson <robj@vmware.com>   
 *
 * ============================================================================
 */

#ifndef _HASHUTIL_H_
#define _HASHUTIL_H_

#include <sys/types.h>
#include <stdlib.h>
#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

// MurmurHash2
uint32_t MurmurHash(const void *buf, size_t length, uint32_t seed);
uint64_t MurmurHash64B ( const void * key, int len, unsigned int
 												seed );
uint64_t MurmurHash64A ( const void * key, int len, unsigned int
 												seed );
__uint128_t MurmurHash128A ( const void * key, int len, unsigned
																		int seed, unsigned int seed2 );

// AES hash
uint64_t AES_HASH(uint64_t x);

uint64_t hash_64(uint64_t key, uint64_t mask);
uint64_t hash_64i(uint64_t key, uint64_t mask);

#ifdef __cplusplus
}
#endif

#endif  // #ifndef _HASHUTIL_H_


