/*
 * ============================================================================
 *
 *       Filename:  cqf.h
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  2017-10-26 11:50:04 AM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Prashant Pandey (), ppandey@cs.stonybrook.edu
 *   Organization:  Stony Brook University
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
#include <aio.h>

#include "cqf/gqf.h"

#define NUM_HASH_BITS 40
#define NUM_Q_BITS 34
#define PAGE_DROP_GRANULARITY (1ULL << 21)
#define BUFFER_SIZE (1ULL << 21)

template <class key_obj>
class CQF {
	public:
		CQF();
		CQF(uint32_t seed);
		CQF(std::string& filename, bool flag);
		CQF(const CQF<key_obj>& copy_cqf);

		void insert(const key_obj& k);

		/* Will return the count. */
		uint64_t query(const key_obj& k);

		void serialize(std::string filename) {
			qf_serialize(&cqf, filename.c_str());
		}

		uint64_t range(void) const { return cqf.metadata->range; }
		uint32_t seed(void) const { return cqf.metadata->seed; }
		uint64_t size(void) const { return cqf.metadata->ndistinct_elts; }
		//uint64_t set_size(void) const { return set.size(); }
		void reset(void) { qf_reset(&cqf); }

		void dump_metadata(void) const { DEBUG_DUMP(&cqf); }

		void drop_pages(uint64_t cur);

		class Iterator {
			public:
				Iterator(QFi it, uint32_t cutoff);
				key_obj operator*(void) const;
				void operator++(void);
				bool done(void) const;

			private:
				/* global buffer to perform read ahead */
				unsigned char buffer[BUFFER_SIZE];
				QFi iter;
				uint32_t cutoff;
				struct aiocb aiocb;
		};

		Iterator begin(uint32_t cutoff) const;
		Iterator end(uint32_t cutoff) const;

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

		uint64_t key;
		uint64_t value;
		uint64_t count;
};

template <class key_obj>
CQF<key_obj>::CQF() {
	qf_init(&cqf, 1ULL << 6, NUM_HASH_BITS, 0, true, "", 23423);
}

template <class key_obj>
CQF<key_obj>::CQF(uint32_t seed) {
	qf_init(&cqf, 1ULL << NUM_Q_BITS, NUM_HASH_BITS, 0, true, "", seed);
}

template <class key_obj>
CQF<key_obj>::CQF(std::string& filename, bool flag) {
	if (flag)
		qf_read(&cqf, filename.c_str());
	else
		qf_deserialize(&cqf, filename.c_str());
}

template <class key_obj>
CQF<key_obj>::CQF(const CQF<key_obj>& copy_cqf) {
	memcpy(cqf, copy_cqf.get_cqf(), sizeof(QF));
}

template <class key_obj>
void CQF<key_obj>::insert(const key_obj& k) {
	qf_insert(&cqf, k.key, k.value, k.count, NO_LOCK);
	// To validate the CQF
	//set.insert(k.key);
}

template <class key_obj>
uint64_t CQF<key_obj>::query(const key_obj& k) {
	return qf_count_key_value(&cqf, k.key, k.value);
}

template <class key_obj>
CQF<key_obj>::Iterator::Iterator(QFi it, uint32_t cutoff)
	: iter(it), cutoff(cutoff) {};

template <class key_obj>
key_obj CQF<key_obj>::Iterator::operator*(void) const {
	uint64_t key = 0, value = 0, count = 0;
	qfi_get(&iter, &key, &value, &count);
	return key_obj(key, value, count);
}

template<class key_obj>
void CQF<key_obj>::Iterator::operator++(void) {
	// Read next 2M bytes from the file offset.
	if (iter.current % PAGE_DROP_GRANULARITY == 0 &&
			iter.current < iter.qf->metadata->nslots - 2 * BUFFER_SIZE) {
		memset(&aiocb, 0, sizeof(struct aiocb));
		aiocb.aio_fildes = iter.qf->mem->fd;
		aiocb.aio_buf = (volatile void*)buffer;
		aiocb.aio_nbytes = BUFFER_SIZE;
		aiocb.aio_offset = (__off_t)(qf_get_addr(iter.qf, iter.current) + BUFFER_SIZE);
		uint32_t ret = aio_read(&aiocb);
		if (ret) {
			std::cerr << "aio_read failed at " << iter.current << " total size " <<
				iter.qf->metadata->nslots << std::endl;
			perror("aio_read");
		}
	}

	qfi_next(&iter);
	// Skip past the cutoff
	do {
		uint64_t key = 0, value = 0, count = 0;
		qfi_get(&iter, &key, &value, &count);
		if (count < cutoff)
			qfi_next(&iter);
		else
			break;
	} while(!qfi_end(&iter));
	// drop pages of the last million slots.
	//static uint64_t last_marker = 1;
	//if (iter.current / PAGE_DROP_GRANULARITY > last_marker + 1) {
		//uint64_t start_idx = last_marker * PAGE_DROP_GRANULARITY;
		//uint64_t end_idx = (last_marker + 1) * PAGE_DROP_GRANULARITY;
		//qf_drop_pages(iter.qf, start_idx, end_idx);
		//last_marker += 1;
	//}
}

template<class key_obj>
bool CQF<key_obj>::Iterator::done(void) const {
	return qfi_end(&iter);
}

template<class key_obj>
typename CQF<key_obj>::Iterator CQF<key_obj>::begin(uint32_t cutoff) const {
	QFi qfi;
	qf_iterator(&this->cqf, &qfi, 0);
	// Skip past the cutoff
	do {
		uint64_t key = 0, value = 0, count = 0;
		qfi_get(&qfi, &key, &value, &count);
		if (count < cutoff)
			qfi_next(&qfi);
		else
			break;
	} while(!qfi_end(&qfi));

	return Iterator(qfi, cutoff);
}

template<class key_obj>
typename CQF<key_obj>::Iterator CQF<key_obj>::end(uint32_t cutoff) const {
	QFi qfi;
	qf_iterator(&this->cqf, &qfi, 0xffffffffffffffff);
	return Iterator(qfi, cutoff);
}

#endif
