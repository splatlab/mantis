/*
 * ============================================================================
 *
 *        Authors:  Prashant Pandey <ppandey@cs.stonybrook.edu>
 *                  Rob Johnson <robj@vmware.com>
 *                  Rob Patro (rob.patro@cs.stonybrook.edu)
 *
 * ============================================================================
 */

#ifndef _CQF_H_
#define _CQF_H_

#include <iostream>
#include <cassert>
#include <unordered_set>

#include <inttypes.h>
#include <string.h>
#include <math.h>
#include <sys/mman.h>

#include "gqf/gqf.h"
#include "gqf/gqf_int.h"
#include "gqf/gqf_file.h"
#include "util.h"

#define NUM_HASH_BITS 14
#define NUM_Q_BITS 6
#define SEED 2038074761

enum readmode {
	CQF_MMAP,
	CQF_FREAD
};

template <class key_obj>
class CQF {
	public:
		CQF();
		CQF(uint64_t q_bits, uint64_t key_bits, enum qf_hashmode hash, uint32_t
				seed);
		CQF(uint64_t q_bits, uint64_t key_bits, enum qf_hashmode hash, uint32_t
				seed, std::string filename);
		CQF(std::string& filename, enum readmode flag);
		CQF(const CQF<key_obj>& copy_cqf);

		int insert(const key_obj& k, uint8_t flags);

		/* Will return the count. */
		uint64_t query(const key_obj& k, uint8_t flags);

		uint64_t inner_prod(const CQF<key_obj>& in_cqf);

		void serialize(std::string filename) {
			qf_serialize(&cqf, filename.c_str());
		}

		void close() { qf_closefile(&cqf); }
		void delete_file() { qf_deletefile(&cqf); }

		void set_auto_resize(void) { qf_set_auto_resize(&cqf, true); }
		int64_t get_unique_index(const key_obj& k, uint8_t flags) const {
			return qf_get_unique_index(&cqf, k.key, k.value, flags);
		}

		bool is_exact(void) const;
		enum qf_hashmode hash_mode(void) const { return cqf.metadata->hash_mode; }
		bool check_similarity(const CQF *other_cqf) const;

		const QF* get_cqf(void) const { return &cqf; }
		uint64_t range(void) const { return cqf.metadata->range; }
		uint32_t seed(void) const { return cqf.metadata->seed; }
		uint64_t numslots(void) const { return cqf.metadata->nslots; }
		uint32_t keybits(void) const { return cqf.metadata->key_bits; }
		uint64_t total_elts(void) const { return cqf.metadata->nelts; }
		uint64_t dist_elts(void) const { return cqf.metadata->ndistinct_elts; }
		//uint64_t set_size(void) const { return set.size(); }
		void reset(void) { qf_reset(&cqf); }

		void dump_metadata(void) const { qf_dump_metadata(&cqf); }

		void drop_pages(uint64_t cur);

		class Iterator {
			public:
				Iterator(QFi it);
				Iterator(QFi it, uint64_t endHash);
				key_obj operator*(void) const;
				void operator++(void);
				bool done(void) const;
				bool reachedHashLimit(void) const;

				key_obj get_cur_hash(void) const;

				QFi iter;
			private:
				uint64_t endHash;
		};

		Iterator begin(void) const;
		Iterator end(void) const;
		Iterator setIteratorLimits(uint64_t start_hash, uint64_t end_hash) const;

	private:
		QF cqf;
		//std::unordered_set<uint64_t> set;
};

class KeyObject {
	public:
		KeyObject() : key(0), value(0), count(0) {};

		KeyObject(uint64_t k, uint64_t v, uint64_t c) : key(k),
		value(v), count(c) {};

		KeyObject(const KeyObject& k) : key(k.key), value(k.value), count(k.count) {};

		bool operator==(KeyObject k) { return key == k.key; }

		typedef uint64_t kmer_t;
		kmer_t key;
		uint64_t value;
		uint64_t count;
};

template <class key_obj>
CQF<key_obj>::CQF() {
	if (!qf_malloc(&cqf, 1ULL << NUM_Q_BITS, NUM_HASH_BITS, 0, QF_HASH_DEFAULT,
								 SEED)) {
		ERROR("Can't allocate the CQF");
		exit(EXIT_FAILURE);
	}
}

template <class key_obj>
CQF<key_obj>::CQF(uint64_t q_bits, uint64_t key_bits, enum qf_hashmode hash,
									uint32_t seed) {
	if (!qf_malloc(&cqf, 1ULL << q_bits, key_bits, 0, hash, seed)) {
		ERROR("Can't allocate the CQF");
		exit(EXIT_FAILURE);
	}
}

template <class key_obj>
CQF<key_obj>::CQF(uint64_t q_bits, uint64_t key_bits, enum qf_hashmode hash,
									uint32_t seed, std::string filename) {
	if (!qf_initfile(&cqf, 1ULL << q_bits, key_bits, 0, hash, seed,
									 filename.c_str())) {
		ERROR("Can't allocate the CQF");
		exit(EXIT_FAILURE);
	}
}

template <class key_obj>
CQF<key_obj>::CQF(std::string& filename, enum readmode flag) {
	uint64_t size = 0;
	if (flag == CQF_MMAP)
	 size = qf_usefile(&cqf, filename.c_str(), QF_USEFILE_READ_ONLY);
	else if (flag == CQF_FREAD)
		size = qf_deserialize(&cqf, filename.c_str());
	else {
		ERROR("Wrong CQF read mode.");
		exit(EXIT_FAILURE);
	}

	if (size == 0) {
		ERROR("Can't read/deserialize the CQF");
		exit(EXIT_FAILURE);
	}
}

template <class key_obj> CQF<key_obj>::CQF(const CQF<key_obj>& copy_cqf) {
	memcpy(reinterpret_cast<void*>(&cqf),
				 reinterpret_cast<void*>(const_cast<QF*>(&copy_cqf.cqf)), sizeof(QF));
}

template <class key_obj>
int CQF<key_obj>::insert(const key_obj& k, uint8_t flags) {
	return qf_insert(&cqf, k.key, k.value, k.count, flags);
	// To validate the CQF
	//set.insert(k.key);
}

template <class key_obj>
uint64_t CQF<key_obj>::query(const key_obj& k, uint8_t flags) {
	return qf_count_key_value(&cqf, k.key, k.value, flags);
}

template <class key_obj>
uint64_t CQF<key_obj>::inner_prod(const CQF<key_obj>& in_cqf) {
	return qf_inner_product(&cqf, in_cqf.get_cqf());
}

template <class key_obj>
bool CQF<key_obj>::is_exact(void) const {
	if (cqf.metadata->hash_mode == QF_HASH_INVERTIBLE)
		return true;
	return false;
}

template <class key_obj>
bool CQF<key_obj>::check_similarity(const CQF *other_cqf) const {
	if (hash_mode() != other_cqf->hash_mode() || seed() != other_cqf->seed() ||
			keybits() != other_cqf->keybits() || range() != other_cqf->range())
		return false;
	return true;
}

template <class key_obj>
CQF<key_obj>::Iterator::Iterator(QFi it)
	: iter(it) {};

template <class key_obj>
CQF<key_obj>::Iterator::Iterator(QFi it, uint64_t endHashIn)
		: iter(it), endHash(endHashIn) {};

template <class key_obj>
key_obj CQF<key_obj>::Iterator::operator*(void) const {
	uint64_t key = 0, value = 0, count = 0;
	qfi_get_key(&iter, &key, &value, &count);
	return key_obj(key, value, count);
}

template <class key_obj>
key_obj CQF<key_obj>::Iterator::get_cur_hash(void) const {
	uint64_t key = 0, value = 0, count = 0;
	qfi_get_hash(&iter, &key, &value, &count);
	return key_obj(key, value, count);
}

template<class key_obj>
void CQF<key_obj>::Iterator::operator++(void) {
	qfi_next(&iter);
}

/* Currently, the iterator only traverses forward. So, we only need to check
 * the right side limit.
 */
template<class key_obj>
bool CQF<key_obj>::Iterator::done(void) const {
	return qfi_end(&iter);
}

/* Currently, the iterator only traverses forward. So, we only need to check
 * the right side limit.
 */
template<class key_obj>
bool CQF<key_obj>::Iterator::reachedHashLimit(void) const {
	uint64_t key = 0, value = 0, count = 0;
	qfi_get_hash(&iter, &key, &value, &count);
	return key >= endHash || qfi_end(&iter);
}

template<class key_obj>
typename CQF<key_obj>::Iterator CQF<key_obj>::begin(void) const {
	QFi qfi;
	qf_iterator_from_position(&this->cqf, &qfi, 0);
	return Iterator(qfi);
}

template<class key_obj>
typename CQF<key_obj>::Iterator CQF<key_obj>::end(void) const {
	QFi qfi;
	qf_iterator_from_position(&this->cqf, &qfi, 0xffffffffffffffff);
	return Iterator(qfi, UINT64_MAX);
}

template<class key_obj>
typename CQF<key_obj>::Iterator CQF<key_obj>::setIteratorLimits(uint64_t start_hash, uint64_t end_hash) const {
	QFi qfi;
	qf_iterator_from_key_value(&this->cqf, &qfi, start_hash, 0, QF_KEY_IS_HASH);
	return Iterator(qfi, end_hash);
}


#endif
