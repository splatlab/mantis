#include "bitvector.h"

BitVector::BitVector(uint64_t size) : size(size) {
		bits = sdsl::bit_vector(size);
}

void BitVector::reset() {
	bits = sdsl::bit_vector(size);
}

bool BitVector::operator[](uint64_t idx) {
	assert(idx < size);
	return bits[idx];
}

void BitVector::set(uint64_t idx) {
	assert(idx < size);
	bits[idx] = 1;
}

void BitVector::resize(const uint64_t len) {
	bits.bit_resize(len);
	size = len;
}

BitVectorRRR::BitVectorRRR(std::string& filename) {
	sdsl::load_from_file(rrr_bits, filename);
	size = rrr_bits.size();
	DEBUG_CDBG("Read rrr bit vector of size " << size << " from file " <<
						 filename);
}

bool BitVectorRRR::operator[](uint64_t idx) {
	assert(idx < size);
	return rrr_bits[idx];
}

bool BitVectorRRR::serialize(std::string& filename) {
	DEBUG_CDBG("Serializing rrr bit vector of size " << size << " to file " <<
						 filename);
	return sdsl::store_to_file(rrr_bits, filename);
}

