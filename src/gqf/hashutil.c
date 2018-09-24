/* -*- Mode: C++; c-basic-offset: 4; indent-tabs-mode: nil -*- */
// Pulled from lookup3.c by Bob Jenkins
#include "gqf/hashutil.h"

//-----------------------------------------------------------------------------
// MurmurHash2, by Austin Appleby
// Note - This code makes a few assumptions about how your machine behaves -
// 1. We can read a 4-byte value from any address without crashing
// 2. sizeof(int) == 4
// And it has a few limitations -
// 1. It will not work incrementally.
// 2. It will not produce the same results on little-endian and big-endian
//    machines.
// All code is released to the public domain. For business purposes,
// Murmurhash is under the MIT license.


uint32_t MurmurHash(const void* buf, size_t len, uint32_t seed)
{
	// 'm' and 'r' are mixing constants generated offline.
	// They're not really 'magic', they just happen to work well.

	const unsigned int m = 0x5bd1e995;
	const int r = 24;

	// Initialize the hash to a 'random' value
	uint32_t h = seed ^ len;

	// Mix 4 bytes at a time into the hash
	const unsigned char * data = (const unsigned char *)buf;

	while(len >= 4) {
		unsigned int k = *(unsigned int *)data;

		k *= m;
		k ^= k >> r;
		k *= m;

		h *= m;
		h ^= k;

		data += 4;
		len -= 4;
	}

	// Handle the last few bytes of the input array
	switch(len) {
		case 3: h ^= data[2] << 16;
		case 2: h ^= data[1] << 8;
		case 1: h ^= data[0];
						h *= m;
	};

	// Do a few final mixes of the hash to ensure the last few
	// bytes are well-incorporated.
	h ^= h >> 13;
	h *= m;
	h ^= h >> 15;
	return h;
}

// Thomas Wang's integer hash functions. See
// <https://gist.github.com/lh3/59882d6b96166dfc3d8d> for a snapshot.

uint64_t hash_64(uint64_t key, uint64_t mask)
{
	key = (~key + (key << 21)) & mask; // key = (key << 21) - key - 1;
	key = key ^ key >> 24;
	key = ((key + (key << 3)) + (key << 8)) & mask; // key * 265
	key = key ^ key >> 14;
	key = ((key + (key << 2)) + (key << 4)) & mask; // key * 21
	key = key ^ key >> 28;
	key = (key + (key << 31)) & mask;
	return key;
}

// The inversion of hash_64(). Modified from
// <https://naml.us/blog/tag/invertible>
uint64_t hash_64i(uint64_t key, uint64_t mask)
{
	uint64_t tmp;

	// Invert key = key + (key << 31)
	tmp = (key - (key << 31));
	key = (key - (tmp << 31)) & mask;

	// Invert key = key ^ (key >> 28)
	tmp = key ^ key >> 28;
	key = key ^ tmp >> 28;

	// Invert key *= 21
	key = (key * 14933078535860113213ull) & mask;

	// Invert key = key ^ (key >> 14)
	tmp = key ^ key >> 14;
	tmp = key ^ tmp >> 14;
	tmp = key ^ tmp >> 14;
	key = key ^ tmp >> 14;

	// Invert key *= 265
	key = (key * 15244667743933553977ull) & mask;

	// Invert key = key ^ (key >> 24)
	tmp = key ^ key >> 24;
	key = key ^ tmp >> 24;

	// Invert key = (~key) + (key << 21)
	tmp = ~key;
	tmp = ~(key - (tmp << 21));
	tmp = ~(key - (tmp << 21));
	key = ~(key - (tmp << 21)) & mask;

	return key;
}

__uint128_t MurmurHash128A ( const void * key, int len,
																			 unsigned int seed1, unsigned int
																			 seed2 ) {
	__uint128_t ret_hash;
	ret_hash = MurmurHash64A(key, len, seed1);
	ret_hash = ret_hash << 64;
	ret_hash = ret_hash | MurmurHash64A(key, len, seed2);

	return ret_hash;
}

//-----------------------------------------------------------------------------
// MurmurHash2, 64-bit versions, by Austin Appleby

// The same caveats as 32-bit MurmurHash2 apply here - beware of alignment 
// and endian-ness issues if used across multiple platforms.


// 64-bit hash for 64-bit platforms

uint64_t MurmurHash64A ( const void * key, int len, unsigned int seed )
{
	const uint64_t m = 0xc6a4a7935bd1e995;
	const int r = 47;

	uint64_t h = seed ^ (len * m);

	const uint64_t * data = (const uint64_t *)key;
	const uint64_t * end = data + (len/8);

	while(data != end)
	{
		uint64_t k = *data++;

		k *= m; 
		k ^= k >> r; 
		k *= m; 

		h ^= k;
		h *= m; 
	}

	const unsigned char * data2 = (const unsigned char*)data;

	switch(len & 7)
	{
		case 7: h ^= (uint64_t)data2[6] << 48;
		case 6: h ^= (uint64_t)data2[5] << 40;
		case 5: h ^= (uint64_t)data2[4] << 32;
		case 4: h ^= (uint64_t)data2[3] << 24;
		case 3: h ^= (uint64_t)data2[2] << 16;
		case 2: h ^= (uint64_t)data2[1] << 8;
		case 1: h ^= (uint64_t)data2[0];
						h *= m;
	};

	h ^= h >> r;
	h *= m;
	h ^= h >> r;

	return h;
}


// 64-bit hash for 32-bit platforms

uint64_t MurmurHash64B ( const void * key, int len, unsigned int seed )
{
	const unsigned int m = 0x5bd1e995;
	const int r = 24;

	unsigned int h1 = seed ^ len;
	unsigned int h2 = 0;

	const unsigned int * data = (const unsigned int *)key;

	while(len >= 8)
	{
		unsigned int k1 = *data++;
		k1 *= m; k1 ^= k1 >> r; k1 *= m;
		h1 *= m; h1 ^= k1;
		len -= 4;

		unsigned int k2 = *data++;
		k2 *= m; k2 ^= k2 >> r; k2 *= m;
		h2 *= m; h2 ^= k2;
		len -= 4;
	}

	if(len >= 4)
	{
		unsigned int k1 = *data++;
		k1 *= m; k1 ^= k1 >> r; k1 *= m;
		h1 *= m; h1 ^= k1;
		len -= 4;
	}

	switch(len)
	{
		case 3: h2 ^= ((unsigned char*)data)[2] << 16;
		case 2: h2 ^= ((unsigned char*)data)[1] << 8;
		case 1: h2 ^= ((unsigned char*)data)[0];
						h2 *= m;
	};

	h1 ^= h2 >> 18; h1 *= m;
	h2 ^= h1 >> 22; h2 *= m;
	h1 ^= h2 >> 17; h1 *= m;
	h2 ^= h1 >> 19; h2 *= m;

	uint64_t h = h1;

	h = (h << 32) | h2;

	return h;
}

uint64_t AES_HASH(uint64_t x)
{
	const uint64_t round_keys[32] =
	{ // These were generated by hashing some randomly chosen files on my laptop
		0x795e15dc8136095f, 0x562371660e56b023,
		0x086bb301d2fb5e87, 0x1fe74f801c68d829,
		0x38a19379fd013357, 0x4a7ef2fca0f840f5,
		0x7d2a08bc58553aef, 0x092cfe1997ab8b53,
		0xd18a0c07dac143d4, 0x64e345ef125a576c,
		0x82807902d8211a1f, 0x6985dc4ddcdaf85d,
		0x2214ff750cf750af, 0xb574b4138eb8a37e,
		0x83e11205e8050dd5, 0x2d62b24118df61eb,
		0x8a16453f8f6b6fa1, 0x260c9e8491474d4f,
		0x06eb44d6042ca8ae, 0x43efbd457306b135,
		0xbfcb7ac89f346686, 0xd00362f30651d0d0,
		0x016d3080768968d5, 0x74b4c2e46ef801de,
		0xf623864a4396fe74, 0x9fc26ea69dad6067,
		0xd0eb2f4e08564d99, 0x408b357725ae0297,
		0xd19efb8e82d22151, 0x58c5ead61b7ecc15,
		0x14e904bc8de1c705, 0x1ef79cd4f487912d
	};
	__uint128_t *rks = (__uint128_t *)round_keys;
	uint64_t output;

	asm("movq       %[input],       %%xmm15;"
			"pxor       %[round_keys0], %%xmm15;"
			"aesenc     %[round_keys1], %%xmm15;"
			"aesenc     %[round_keys2], %%xmm15;"
			"aesenc     %[round_keys3], %%xmm15;"
			"aesenc     %[round_keys4], %%xmm15;"
			"aesenc     %[round_keys5], %%xmm15;"
			"aesenc     %[round_keys6], %%xmm15;"
			"aesenc     %[round_keys7], %%xmm15;"
			"aesenc     %[round_keys8], %%xmm15;"
			"aesenc     %[round_keys9], %%xmm15;"
			"aesenclast %[round_keysa], %%xmm15;"
			"vmovq      %%xmm15,        %[output]"
			: [output] "=irm" (output)
			: [input] "irm" (x),
			[round_keys0] "m" (rks[0]),
			[round_keys1] "m" (rks[1]),
			[round_keys2] "m" (rks[2]),
			[round_keys3] "m" (rks[3]),
			[round_keys4] "m" (rks[4]),
			[round_keys5] "m" (rks[5]),
			[round_keys6] "m" (rks[6]),
			[round_keys7] "m" (rks[7]),
			[round_keys8] "m" (rks[8]),
			[round_keys9] "m" (rks[9]),
			[round_keysa] "m" (rks[10])
				 : "xmm15"
					 );

	return output;
}

