/*
 * ============================================================================
 *
 *       Filename:  coloreddbg.h
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  2017-10-24 08:49:22 PM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Prashant Pandey (), ppandey@cs.stonybrook.edu
 *   Organization:  Stony Brook University
 *
 * ============================================================================
 */

#ifndef _COLORED_DBG_H_
#define _COLORED_DBG_H_

#include <iostream>
#include <fstream>
#include <unordered_map>
#include <set>
#include <unordered_set>

#include <inttypes.h>

#include "libcuckoo/cuckoohash_map.hh"
#include "sparsepp/spp.h"
#include "tsl/sparse_map.h"
#include "sdsl/bit_vectors.hpp"
#include "bitvector.h"
#include "cqf.h"
#include "hashutil.h"
#include "common_types.h"

#define INITIAL_EQ_CLASSES 10000
#define CQF_FILE "dbg_cqf.ser"
#define EQCLASS_FILE "eqclass_rrr.cls"
#define SAMPLEID_FILE "sampleid.lst"

template <typename Key, typename Value, typename Hasher>
  using cdbg_bv_map_t = cuckoohash_map<Key, Value, Hasher>;

template <class qf_obj, class key_obj>
class ColoredDbg {
  	public:
    ColoredDbg(std::string& cqf_file, std::string& eqclass_file, std::string&
							 sample_file);

		ColoredDbg(uint64_t key_bits, uint32_t seed, uint32_t nqf);
		
		void build_sampleid_map(qf_obj *incqfs);

		cdbg_bv_map_t<BitVector, std::pair<uint64_t,uint64_t>,
		sdslhash<BitVector>>& construct(qf_obj *incqfs,
																		std::unordered_map<std::string, uint64_t>&
																		cutoffs, cdbg_bv_map_t<BitVector,
																		std::pair<uint64_t, uint64_t>,
																		sdslhash<BitVector>>& map, uint64_t
																		start_hash, uint64_t end_hash, uint64_t
																		num_kmers);

		const CQF<key_obj> *get_cqf(void) const { return &dbg; }
		const BitVectorRRR get_bitvector(void) const { return eqclasses; }
		uint64_t get_num_eqclasses(void) const { return num_eq_classes.load(); }
		uint64_t get_num_samples(void) const { return num_samples; }
		std::string get_sample(uint32_t id) const;
		uint32_t seed(void) const { return dbg.seed(); }
		uint64_t range(void) const { return dbg.range(); }

		std::unordered_map<uint64_t, uint64_t>
			find_samples(const mantis::QuerySet& kmers);

		void serialize(std::string prefix);

	private:
		void increment_num_eqclasses() { ++num_eq_classes; }
		void add_kmer(key_obj& hash, BitVector& vector);
		void add_eq_class(BitVector vector, uint64_t id);
		uint64_t get_next_available_id(void);

		std::unordered_map<uint64_t, std::string> sampleid_map;
		// bit_vector --> <eq_class_id, abundance>
		cdbg_bv_map_t<BitVector, std::pair<uint64_t, uint64_t>,
			sdslhash<BitVector>> eqclass_map;
		CQF<key_obj> dbg;
		BitVectorRRR eqclasses;
		uint32_t num_samples;
		std::atomic<int> num_eq_classes;
};

template <class T>
class SampleObject {
	public:
		SampleObject() : obj(), sample_id(), id(0) {};
		SampleObject(T o, std::string& s, uint32_t id) : obj(o),
		sample_id(s), id(id) {};
		SampleObject(const SampleObject& o) : obj(o.obj),
		sample_id(o.sample_id), id(o.id) {} ;

		T obj;
		std::string sample_id;
		uint32_t id;
};

template <class T>
struct compare {
	bool operator()(const SampleObject<T>& lhs, const SampleObject<T>& rhs)
	{
		return lhs.obj.key > rhs.obj.key;
	}
};

template <class qf_obj, class key_obj>
inline uint64_t ColoredDbg<qf_obj, key_obj>::get_next_available_id(void) {
	return get_num_eqclasses() + 1;
}

template <class qf_obj, class key_obj>
std::string ColoredDbg<qf_obj, key_obj>::get_sample(uint32_t id) const {
	auto it = sampleid_map.find(id);
	if (it == sampleid_map.end())
		return std::string();
	else
		return it->second;
}

template <class qf_obj, class key_obj>
void ColoredDbg<qf_obj, key_obj>::add_kmer(key_obj& k, BitVector&
																					 vector) {
	// A kmer (hash) is seen only once during the merge process.
	// So we insert every kmer in the dbg
	auto updatefn = [](std::pair<uint64_t, uint64_t> &val) { ++val.second; };
	bool new_vector = eqclass_map.upsert(vector, updatefn,
												std::pair<uint64_t,uint64_t>(get_next_available_id(),
																										 1));
	if (new_vector) {
		increment_num_eqclasses();
	}

	uint64_t eq_id = eqclass_map.find(vector).first;

	k.count = eq_id;	// we use the count to store the eqclass ids
	dbg.insert(k);

	static uint64_t last_size = 0;
	if (dbg.size() % 10000000 == 0 &&
			dbg.size() != last_size) {
		last_size = dbg.size();
		std::cout << "Kmers merged: " << dbg.size() << " Num eq classes: " <<
			get_num_eqclasses() << std::endl;
	}
}

template <class qf_obj, class key_obj>
std::unordered_map<uint64_t, uint64_t>
ColoredDbg<qf_obj,key_obj>::find_samples(const mantis::QuerySet& kmers) {
	// Find a list of eq classes and the number of kmers that belong those eq
	// classes.
	std::unordered_map<uint64_t, uint64_t> query_eqclass_map;
	for (auto k : kmers) {
		key_obj key(k, 0, 0);
		uint64_t eqclass = dbg.query(key);
		if (eqclass)
			query_eqclass_map[eqclass] += 1;
	}

	std::unordered_map<uint64_t, uint64_t> sample_map;
	for (auto it = query_eqclass_map.begin(); it != query_eqclass_map.end();
			 ++it) {
		auto eqclass_id = it->first;
		auto count = it->second;
		// counter starts from 1.
		uint64_t start_idx = (eqclass_id - 1) * num_samples;
		for (uint32_t w = 0; w <= num_samples / 64; w++) {
			uint32_t len = std::min((uint32_t)64, num_samples - w * 64);
			uint64_t wrd = eqclasses.get_int(start_idx, len);
			for (uint32_t i = 0, sCntr = w * 64; i < len; i++, sCntr++)
					if ((wrd >> i) & 0x01)
							sample_map[sCntr] += count;
			start_idx += len;
		}
	}
	return sample_map;
}

template <class qf_obj, class key_obj>
void ColoredDbg<qf_obj, key_obj>::serialize(std::string prefix) {
	dbg.serialize(prefix + CQF_FILE);
	BitVector final_bv(get_num_eqclasses() * num_samples);
	auto vector_it = eqclass_map.lock_table();
	std::ofstream tmpfile(prefix + "eqclass_dist.lst");
	for (const auto &it : vector_it) {
		tmpfile << it.second.first << " " << it.second.second << std::endl;
		BitVector vector = it.first;
		// counter starts from 1.
		uint64_t cur_id = it.second.first - 1;
		uint64_t start_idx = cur_id * num_samples;
		for (uint32_t i = 0; i < num_samples; i++, start_idx++)
			if (vector[i])
				final_bv.set(start_idx);
	}
	BitVectorRRR final_com_bv(final_bv);
	std::string bv_file(prefix + EQCLASS_FILE);
	final_com_bv.serialize(bv_file);
	std::ofstream opfile(prefix + SAMPLEID_FILE);
	for (auto sample : sampleid_map)
		opfile << sample.first << " " << sample.second << std::endl;
	
	opfile.close();
	tmpfile.close();
}

template <class qf_obj, class key_obj>
cdbg_bv_map_t<BitVector, std::pair<uint64_t, uint64_t>,
	sdslhash<BitVector>>& ColoredDbg<qf_obj, key_obj>::construct(qf_obj *incqfs,
															std::unordered_map<std::string, uint64_t>& cutoffs,
															cdbg_bv_map_t<BitVector, std::pair<uint64_t, uint64_t>,
															sdslhash<BitVector>>& map, uint64_t start_hash,
															uint64_t end_hash, uint64_t num_kmers) {
																		
	uint32_t nqf = num_samples;
	uint64_t counter = 0;

	eqclass_map = map;
	dbg.reset();

	// merge all input CQFs into the final QF
	typename CQF<key_obj>::Iterator *it_incqfs =
		(typename CQF<key_obj>::Iterator*)calloc(num_samples, sizeof(typename
																																 CQF<key_obj>::Iterator));

	// Initialize all iterators with sample specific cutoffs.
	for (uint32_t i = 0; i < num_samples; i++) {
		auto it = cutoffs.find(incqfs[i].sample_id);
		if (it == cutoffs.end()) {
			std::cerr << "Sample id " <<  incqfs[i].sample_id << " not found in" <<
				" cutoff list." << std::endl;
			abort();
		} else
			it_incqfs[i] = incqfs[i].obj->limits(start_hash, end_hash, it->second);
	}

	std::priority_queue<SampleObject<KeyObject>,
		std::vector<SampleObject<KeyObject>>, compare<KeyObject>> minheap;
	// Insert the first key from each CQF in minheap.
	for (uint32_t i = 0; i < num_samples; i++) {
		if (it_incqfs[i].done())
			continue;
		KeyObject key = *it_incqfs[i];
		SampleObject<KeyObject> obj(key, incqfs[i].sample_id, i);
		minheap.push(obj);
	}

	while (!minheap.empty()) {
		assert(minheap.size() == nqf);
		SampleObject<KeyObject> cur;
		BitVector eq_class(num_samples);
		// Get the smallest key from minheap and update the eqclass vector
		cur = minheap.top();
		eq_class.set(cur.id);
		minheap.pop();
		// Keep poping keys from minheap until you see a different key.
		// While poping keys build the eq class for cur.
		// Increment iterators for all CQFs whose keys are popped.
		while (!minheap.empty() && cur.obj.key == minheap.top().obj.key) {
			uint32_t id = minheap.top().id;
			eq_class.set(id);
			minheap.pop();
			++it_incqfs[id];
			if (it_incqfs[id].done())	// If the iterator is done then decrement nqf
				nqf--;
			else {	// Insert the current iterator head in minHeap
				KeyObject key = *it_incqfs[id];
				SampleObject<KeyObject> obj(key, incqfs[id].sample_id, id);
				minheap.push(obj);
			}
		}
		// Move the iterator of the smallest key.
		++it_incqfs[cur.id];
		if (it_incqfs[cur.id].done())	// If the iterator is done then decrement nqf
			nqf--;
		else {	// Insert the current iterator head in minHeap
			KeyObject key = *it_incqfs[cur.id];
			SampleObject<KeyObject> obj(key, incqfs[cur.id].sample_id, cur.id);
			minheap.push(obj);
		}
		// Add <kmer, vector> in the cdbg
		add_kmer(cur.obj, eq_class);
		counter++;
		if (counter > num_kmers)
			break;
	}

	return eqclass_map;
}

template <class qf_obj, class key_obj>
void ColoredDbg<qf_obj, key_obj>::build_sampleid_map(qf_obj *incqfs) {
	for (uint32_t i = 0; i < num_samples; i++) {
		std::pair<uint32_t, std::string> pair(incqfs[i].id, incqfs[i].sample_id);
		sampleid_map.insert(pair);
	}
}

template <class qf_obj, class key_obj>
ColoredDbg<qf_obj, key_obj>::ColoredDbg(uint64_t key_bits, uint32_t seed,
																				uint32_t nqf) :
	dbg(key_bits, seed), eqclasses(nqf * INITIAL_EQ_CLASSES),
	num_samples(nqf), num_eq_classes(0) {}

template <class qf_obj, class key_obj>
ColoredDbg<qf_obj, key_obj>::ColoredDbg(std::string& cqf_file, std::string&
																				eqclass_file, std::string& sample_file)
: dbg(cqf_file, false), eqclasses(eqclass_file) {
	num_samples = 0;
	std::ifstream sampleid(sample_file.c_str());
	std::string sample;
	uint32_t id;
	while (sampleid >> id >> sample) {
		std::pair<uint32_t, std::string> pair(id, sample);
		sampleid_map.insert(pair);
		num_samples++;
	}
	sampleid.close();
}

#endif
