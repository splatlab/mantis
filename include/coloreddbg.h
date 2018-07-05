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
#include <chrono>

#include <inttypes.h>

#include "sparsepp/spp.h"
#include "tsl/sparse_map.h"
#include "sdsl/bit_vectors.hpp"
#include "bitvector.h"
#include "cqf.h"
#include "hashutil.h"
#include "common_types.h"
#include "mantisconfig.hpp"

struct hash128 {
	uint64_t operator()(const __uint128_t& val128) const
	{
		__uint128_t val = val128;
		// Using the same seed as we use in k-mer hashing.
		return HashUtil::MurmurHash64A((void*)&val, sizeof(__uint128_t),
																	 2038074743);
	}
};

template <typename Key, typename Value>
  using cdbg_bv_map_t = spp::sparse_hash_map<Key, Value, hash128>;

using default_cdbg_bv_map_t = cdbg_bv_map_t<__uint128_t, std::pair<uint64_t,uint64_t>>;

template <class qf_obj, class key_obj>
class ColoredDbg {
  	public:
		ColoredDbg(std::string& cqf_file, std::vector<std::string>& eqclass_files,
							 std::string& sample_file);

		ColoredDbg(uint64_t qbits, uint64_t key_bits, uint32_t seed,
							 std::string& prefix, uint64_t nqf);
		
		void build_sampleid_map(qf_obj *incqfs);

    default_cdbg_bv_map_t&
			construct(qf_obj *incqfs, default_cdbg_bv_map_t& map, uint64_t num_kmers);

		void set_console(spdlog::logger* c) { console = c; }
		const CQF<key_obj> *get_cqf(void) const { return &dbg; }
		uint64_t get_num_bitvectors(void) const;
		uint64_t get_num_eqclasses(void) const { return eqclass_map.size(); }
		uint64_t get_num_samples(void) const { return num_samples; }
		std::string get_sample(uint32_t id) const;
		uint32_t seed(void) const { return dbg.seed(); }
		uint64_t range(void) const { return dbg.range(); }

		std::unordered_map<uint64_t, uint64_t>
			find_samples(const mantis::QuerySet& kmers);

		void serialize();
		void reinit(default_cdbg_bv_map_t& map);
		void set_flush_eqclass_dist(void) { flush_eqclass_dis = true; }

	private:
    // returns true if adding this k-mer increased the number of equivalence classes
    // and false otherwise.
		bool add_kmer(key_obj& hash, BitVector& vector);
		void add_bitvector(BitVector& vector, uint64_t eq_id);
		void add_eq_class(BitVector vector, uint64_t id);
		uint64_t get_next_available_id(void);
		void bv_buffer_serialize();
		void reshuffle_bit_vectors(cdbg_bv_map_t<__uint128_t, std::pair<uint64_t,
															 uint64_t>>& map);

		std::unordered_map<uint64_t, std::string> sampleid_map;
		// bit_vector --> <eq_class_id, abundance>
		cdbg_bv_map_t<__uint128_t, std::pair<uint64_t, uint64_t>> eqclass_map;
		CQF<key_obj> dbg;
		BitVector bv_buffer;
		std::vector<BitVectorRRR> eqclasses;
		std::string prefix;
		uint64_t num_samples;
		uint64_t num_serializations;
		bool flush_eqclass_dis{false};
    std::time_t start_time_;
		spdlog::logger* console;
};

template <class T>
class SampleObject {
	public:
		SampleObject() : obj(), cutoff(0), sample_id(), id(0) {};
		SampleObject(T o, uint32_t c = 0, std::string& s = std::string(),
								 uint32_t id = 0) : obj(o), cutoff(c), sample_id(s), id(id) {};
		SampleObject(const SampleObject& o) : obj(o.obj), cutoff(o.cutoff),
		sample_id(o.sample_id), id(o.id) {} ;

		T obj;
		uint32_t cutoff;
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
uint64_t ColoredDbg<qf_obj, key_obj>::get_num_bitvectors(void) const {
	uint64_t total = 0;
	for (uint32_t i = 0; i < num_serializations; i++)
		total += eqclasses[i].bit_size();

	return total / num_samples;
}

template <class qf_obj, class key_obj>
void ColoredDbg<qf_obj,
		 key_obj>::reshuffle_bit_vectors(cdbg_bv_map_t<__uint128_t,
																		 std::pair<uint64_t, uint64_t>>& map) {
	BitVector new_bv_buffer(mantis::NUM_BV_BUFFER * num_samples);
	for (auto& it_input : map) {
		auto it_local = eqclass_map.find(it_input.first);
		if (it_local == eqclass_map.end()) {
			console->error("Can't find the vector hash during shuffling");
			exit(1);
		} else {
			assert(it_local->second.first <= mantis::NUM_BV_BUFFER && it_input.second.first
						 <= mantis::NUM_BV_BUFFER);
			uint64_t src_idx = ((it_local->second.first - 1) * num_samples);
			uint64_t dest_idx = ((it_input.second.first - 1) * num_samples);
			for (uint32_t i = 0; i < num_samples; i++, src_idx++, dest_idx++)
				if (bv_buffer[src_idx])
					new_bv_buffer.set(dest_idx);
		}
	}
	bv_buffer = new_bv_buffer;
}

template <class qf_obj, class key_obj>
void ColoredDbg<qf_obj, key_obj>::reinit(cdbg_bv_map_t<__uint128_t,
																		 std::pair<uint64_t, uint64_t>>& map) {
	dbg.reset();
	reshuffle_bit_vectors(map);
	// Check if the current bit vector buffer is full and needs to be serialized.
	// This happens when the sampling phase fills up the bv buffer.
	if (get_num_eqclasses() % mantis::NUM_BV_BUFFER == 0) {
		// The bit vector buffer is full.
		console->info("Serializing bit vector with {} eq classes.",
									get_num_eqclasses());
		bv_buffer_serialize();
	}
	eqclass_map = map;
}

template <class qf_obj, class key_obj>
bool ColoredDbg<qf_obj, key_obj>::add_kmer(key_obj& k, BitVector&
																					 vector) {
	// A kmer (hash) is seen only once during the merge process.
	// So we insert every kmer in the dbg
	uint64_t eq_id = 1;
	__uint128_t vec_hash = HashUtil::MurmurHash128A((void*)vector.data(),
																								 vector.capacity(), 2038074743,
																								 2038074751);

	auto it = eqclass_map.find(vec_hash);
  bool added_eq_class{false};
	// Find if the eqclass of the kmer is already there.
	// If it is there then increment the abundance.
	// Else create a new eq class.
	if (it == eqclass_map.end()) {
		// eq class is seen for the first time.
		eq_id = get_next_available_id();
		eqclass_map.emplace(std::piecewise_construct,
															std::forward_as_tuple(vec_hash),
															std::forward_as_tuple(eq_id, 1));
		add_bitvector(vector, eq_id - 1);
    added_eq_class = true;
	} else { // eq class is seen before so increment the abundance.
		eq_id = it->second.first;
    // with standard map
    it->second.second += 1; // update the abundance.
	}

	k.count = eq_id;	// we use the count to store the eqclass ids
	dbg.insert(k);
  return added_eq_class;
}

template <class qf_obj, class key_obj>
void ColoredDbg<qf_obj, key_obj>::add_bitvector(BitVector& vector, uint64_t
																								eq_id) {
	uint64_t start_idx = (eq_id  % mantis::NUM_BV_BUFFER) * num_samples;
	for (uint32_t i = 0; i < num_samples; i++, start_idx++)
		if (vector[i])
			bv_buffer.set(start_idx);
}

template <class qf_obj, class key_obj>
void ColoredDbg<qf_obj, key_obj>::bv_buffer_serialize() {
	BitVector bv_temp(bv_buffer);
	if (get_num_eqclasses() % mantis::NUM_BV_BUFFER > 0) {
		bv_temp.resize((get_num_eqclasses() % mantis::NUM_BV_BUFFER) * num_samples);
	}

	BitVectorRRR final_com_bv(bv_temp);
	std::string bv_file(prefix + std::to_string(num_serializations) + "_" +
                      mantis::EQCLASS_FILE);
	final_com_bv.serialize(bv_file);
	bv_buffer.reset();
	num_serializations++;
}

template <class qf_obj, class key_obj>
void ColoredDbg<qf_obj, key_obj>::serialize() {
	// serialize the CQF
	dbg.serialize(prefix + mantis::CQF_FILE);

	// serialize the bv buffer last time if needed
	if (get_num_eqclasses() % mantis::NUM_BV_BUFFER > 1)
		bv_buffer_serialize();

	//serialize the eq class id map
	std::ofstream opfile(prefix + mantis::SAMPLEID_FILE);
	for (auto sample : sampleid_map)
		opfile << sample.first << " " << sample.second << std::endl;
	opfile.close();

	if (flush_eqclass_dis) {
		// dump eq class abundance dist for further analysis.
		std::ofstream tmpfile(prefix + "eqclass_dist.lst");
		for (auto sample : eqclass_map)
			tmpfile << sample.second.first << " " << sample.second.second << std::endl;
		tmpfile.close();
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
		uint64_t start_idx = (eqclass_id - 1);
		uint64_t bucket_idx = start_idx / mantis::NUM_BV_BUFFER;
		uint64_t bucket_offset = (start_idx % mantis::NUM_BV_BUFFER) * num_samples;
		for (uint32_t w = 0; w <= num_samples / 64; w++) {
			uint64_t len = std::min((uint64_t)64, num_samples - w * 64);
			uint64_t wrd = eqclasses[bucket_idx].get_int(bucket_offset, len);
			for (uint32_t i = 0, sCntr = w * 64; i < len; i++, sCntr++)
					if ((wrd >> i) & 0x01)
							sample_map[sCntr] += count;
			bucket_offset += len;
		}
	}
	return sample_map;
}

template <class qf_obj, class key_obj>
cdbg_bv_map_t<__uint128_t, std::pair<uint64_t, uint64_t>>& ColoredDbg<qf_obj,
	key_obj>::construct(qf_obj *incqfs, cdbg_bv_map_t<__uint128_t,
											std::pair<uint64_t, uint64_t>>& map, uint64_t num_kmers)
{
	uint32_t nqf = 0;
	uint64_t counter = 0;

  bool is_sampling = (num_kmers < std::numeric_limits<uint64_t>::max());

	// merge all input CQFs into the final QF
  std::vector<typename CQF<key_obj>::Iterator> it_incqfs;
  it_incqfs.reserve(num_samples);

	// Initialize all iterators with sample specific cutoffs.
	for (uint32_t i = 0; i < num_samples; i++) {
		it_incqfs.emplace_back(incqfs[i].obj->begin(incqfs[i].cutoff));
  }

	std::priority_queue<SampleObject<KeyObject>,
		std::vector<SampleObject<KeyObject>>, compare<KeyObject>> minheap;
	// Insert the first key from each CQF in minheap.
	for (uint32_t i = 0; i < num_samples; i++) {
		if (it_incqfs[i].done())
			continue;
		KeyObject key = *it_incqfs[i];
		minheap.emplace(key, incqfs[i].cutoff, incqfs[i].sample_id, i);
		nqf++;
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
				minheap.emplace(key, incqfs[id].cutoff, incqfs[id].sample_id, id);
			}
		}
		// Move the iterator of the smallest key.
		++it_incqfs[cur.id];
		if (it_incqfs[cur.id].done())	// If the iterator is done then decrement nqf
			nqf--;
		else {	// Insert the current iterator head in minHeap
			KeyObject key = *it_incqfs[cur.id];
			minheap.emplace(key, incqfs[cur.id].cutoff, incqfs[cur.id].sample_id, cur.id);
		}
		// Add <kmer, vector> in the cdbg
		bool added_eq_class = add_kmer(cur.obj, eq_class);
		counter++;

		// Progress tracker
		static uint64_t last_size = 0;
		if (dbg.size() % 10000000 == 0 &&
				dbg.size() != last_size) {
			last_size = dbg.size();
			console->info("Kmers merged: {}  Num eq classes: {}  Total time: {}",
										dbg.size(), get_num_eqclasses(), time(nullptr) - start_time_);
		}

		// Check if the bit vector buffer is full and needs to be serialized.
		if (added_eq_class and (get_num_eqclasses() % mantis::NUM_BV_BUFFER == 0)) {
			// Check if the process is in the sampling phase.
			if (is_sampling) {
				break;
      } else {
				// The bit vector buffer is full.
				console->info("Serializing bit vector with {} eq classes.",
											get_num_eqclasses());
				bv_buffer_serialize();
			}
		} else if (counter > num_kmers) {
      // Check if the sampling phase is finished based on the number of k-mers.
			break;
    }
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
ColoredDbg<qf_obj, key_obj>::ColoredDbg(uint64_t qbits, uint64_t key_bits,
																				uint32_t seed, std::string& prefix,
																				uint64_t nqf) :
	dbg(qbits, key_bits, seed), bv_buffer(mantis::NUM_BV_BUFFER * nqf),
    prefix(prefix), num_samples(nqf), num_serializations(0), start_time_(std::time(nullptr)) {}

template <class qf_obj, class key_obj>
ColoredDbg<qf_obj, key_obj>::ColoredDbg(std::string& cqf_file,
																				std::vector<std::string>&
																				eqclass_files, std::string&
																				sample_file)
    : dbg(cqf_file, false), bv_buffer(), start_time_(std::time(nullptr)) {
	num_samples = 0;
	num_serializations = 0;

	std::map<int, std::string> sorted_files;
	for (std::string file : eqclass_files) {
		int id = std::stoi(first_part(last_part(file, '/'), '_'));
		sorted_files[id] = file;
	}

	for (auto file : sorted_files) {
		eqclasses.push_back(BitVectorRRR(file.second));
		num_serializations++;
	}

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
