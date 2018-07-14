/*
 * ============================================================================
 *
 *       Filename:  bitvector.h
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  2017-10-25 01:59:49 PM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Prashant Pandey (), ppandey@cs.stonybrook.edu
 *   Organization:  Stony Brook University
 *
 * ============================================================================
 */

#ifndef _BIT_VECTOR_H_
#define _BIT_VECTOR_H_

#include <iostream>
#include <cassert>

#include <inttypes.h>

#include "util.h"
#include "hashutil.h"
#include "sdsl/bit_vectors.hpp"

class BitVector {
	public:
		BitVector() : bits(0), size(0) {};
		BitVector(uint64_t size);

		sdsl::bit_vector get_bits() const {
			return bits;
		}
		bool operator[](uint64_t idx);
		void set(const uint64_t idx);
		uint64_t capacity(void) const { return bits.capacity() / 8; }
		uint64_t bit_size(void) const { return bits.bit_size(); }
		const uint64_t *data(void) const { return bits.data(); }
		void resize(const uint64_t len);
		uint64_t get_int(uint64_t startP, uint64_t len=64) {return bits.get_int(startP, len);}
		bool operator==(const BitVector& b) const { return bits == b.bits; }

	private:
		sdsl::bit_vector bits;
		uint64_t size;
};

class BitVectorRRR {
	public:
		BitVectorRRR() : rrr_bits(BitVector().get_bits()), size(0) {};
		BitVectorRRR(const BitVector& bv) : rrr_bits(bv.get_bits()),
				size(bv.bit_size()) {};
		BitVectorRRR(std::string& filename);

		bool operator[](uint64_t idx);
		bool serialize(std::string& filename);
		uint64_t bit_size(void) const { return rrr_bits.size(); }
		uint64_t get_int(uint64_t startP, uint64_t len=64) const {
			return rrr_bits.get_int(startP, len);
		}

	private:
		sdsl::rrr_vector<63> rrr_bits;
		uint64_t size;
};

template <class T>
struct sdslhash {
	uint64_t operator()(const T& vector) const
	{
		// Using the same seed as we use in k-mer hashing.
		return HashUtil::MurmurHash64A((void*)vector.data(), vector.capacity(),
																	 2038074743);
	}
};

#endif

