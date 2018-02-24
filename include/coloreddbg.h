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

#include "sparsepp/spp.h"
#include "tsl/sparse_map.h"
#include "sdsl/bit_vectors.hpp"
#include "bitvector.h"
#include "cqf.h"
#include "hashutil.h"
#include "common_types.h"

#define NUM_BV_BUFFER 20000000
#define INITIAL_EQ_CLASSES 10000
#define CQF_FILE "dbg_cqf.ser"
#define EQCLASS_FILE "eqclass_rrr.cls"
#define SAMPLEID_FILE "sampleid.lst"
#define ACCUM_EQCLS_CNT_FILE "accum_eqcls_cnt.lst"

extern uint64_t start_time;


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

template <class qf_obj, class key_obj>
class ColoredDbg {
  	public:
		ColoredDbg(std::string& cqf_file, std::vector<std::string>& eqclass_files,
							 std::string& sample_file);

		ColoredDbg(uint64_t key_bits, uint32_t seed, std::string& prefix, uint64_t
							 nqf);
		
		void build_sampleid_map(qf_obj *incqfs);

		cdbg_bv_map_t<__uint128_t, std::pair<uint64_t,uint64_t>>&
			construct(qf_obj *incqfs, std::unordered_map<std::string, uint64_t>&
								cutoffs, cdbg_bv_map_t<__uint128_t, std::pair<uint64_t,
								uint64_t>>& map, uint64_t num_kmers);

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
		void reinit(cdbg_bv_map_t<__uint128_t, std::pair<uint64_t, uint64_t>>&
								map);

	private:
		void add_kmer(key_obj& hash, BitVector& vector);
		void add_bitvector(BitVector& vector, uint64_t eq_id);
		void add_eq_class(BitVector vector, uint64_t id);
		uint64_t get_next_available_id(void);
		void bv_buffer_serialize();
		void reshuffle_bit_vectors(cdbg_bv_map_t<__uint128_t, std::pair<uint64_t,
															 uint64_t>>& map);
		void copy(BitVector& src, BitVector& dst, uint64_t fromStartIdx, uint64_t toStartIdx);
		std::unordered_map<uint64_t, std::string> sampleid_map;
		// bit_vector --> <eq_class_id, abundance>
		cdbg_bv_map_t<__uint128_t, std::pair<uint64_t, uint64_t>> eqclass_map;
		CQF<key_obj> dbg;
		BitVector bv_buffer;
		std::vector<BitVectorRRR> eqclasses;
		std::string prefix;
		uint64_t num_samples;
		uint64_t num_serializations;
		bool is_sampling{false};
		std::vector<uint64_t> accumEqClsCnt;

		uint64_t eqt{0};
		uint64_t cqft{0};
		uint64_t eq_exportt{0};
		uint64_t cqf_exportt{0};
		uint64_t reshufflet{0};
		uint64_t map_hasht{0};
		uint64_t rrrt{0};
		uint64_t cqf_heapt{0};
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
	uint64_t s = time(NULL);
	BitVector new_bv_buffer(NUM_BV_BUFFER * num_samples);
	for (auto& it_input : map) {
		auto it_local = eqclass_map.find(it_input.first);
		if (it_local == eqclass_map.end()) {
			std::cerr << "Can't find the vector hash during shuffling." << std::endl;
			exit(1);
		} else {
			uint64_t src_idx = ((it_local->second.first - 1) * num_samples) %
				NUM_BV_BUFFER;
			uint64_t dest_idx = ((it_input.second.first - 1) * num_samples) %
				NUM_BV_BUFFER;
			copy(bv_buffer, new_bv_buffer, src_idx, dest_idx);
			/* for (uint32_t i = 0; i < num_samples; i++, src_idx++, dest_idx++)
				if (bv_buffer[src_idx])
					new_bv_buffer.set(dest_idx); */
		}
	}
	bv_buffer = new_bv_buffer;
	reshufflet += time(NULL) - s;
}

template <class qf_obj, class key_obj>
void ColoredDbg<qf_obj, key_obj>::reinit(cdbg_bv_map_t<__uint128_t,
																		 std::pair<uint64_t, uint64_t>>& map) {
	dbg.reset();
	reshuffle_bit_vectors(map);
	uint64_t s = time(NULL);
	eqclass_map = map;
	map_hasht += time(NULL) - s;
}

template <class qf_obj, class key_obj>
void ColoredDbg<qf_obj, key_obj>::add_kmer(key_obj& k, BitVector&
																					 vector) {
	// A kmer (hash) is seen only once during the merge process.
	// So we insert every kmer in the dbg
	uint64_t s = time(NULL);
	
	uint64_t eq_id = 1;
	__uint128_t vec_hash = HashUtil::MurmurHash128A((void*)vector.data(),
																								 vector.capacity(), 2038074743,
																								 2038074751);

	auto it = eqclass_map.find(vec_hash);
	map_hasht += time(NULL) - s;
	// Find if the eqclass of the kmer is already there.
	// If it is there then increment the abundance.
	// Else create a new eq class.
	if (it == eqclass_map.end()) {
		// eq class is seen for the first time.
		eq_id = get_next_available_id();
		s = time(NULL);
		eqclass_map.emplace(std::piecewise_construct,
															std::forward_as_tuple(vec_hash),
															std::forward_as_tuple(eq_id, 1));
		map_hasht += time(NULL) - s;
		
		add_bitvector(vector, eq_id - 1);
		
	} else { // eq class is seen before so increment the abundance.
		eq_id = it->second.first;
    // with standard map
    it->second.second += 1; // update the abundance.
	}

	k.count = eq_id;	// we use the count to store the eqclass ids
	s= time(NULL);
	dbg.insert(k);
	cqft += time(NULL) - s;

	// Serialize bit vectors if buffer is full.
	if (it == eqclass_map.end() && get_num_eqclasses() % NUM_BV_BUFFER == 0) {
		PRINT_CDBG("Serializing bit vector with " << get_num_eqclasses() <<
							 " eq classes.");
		bv_buffer_serialize();
		PRINT_CDBG("Export time: " << (time(NULL) - s));
	}

	s = time(NULL);
	static uint64_t last_size = 0;
	if (dbg.size() % 10000000 == 0 &&
			dbg.size() != last_size) {
		last_size = dbg.size();
		PRINT_CDBG("Kmers merged: " << dbg.size() << " Num eq classes: " <<
			get_num_eqclasses() <<  " Total time: " << time(NULL) - start_time);
	}
	cqft += time(NULL) - s;
}

template <class qf_obj, class key_obj>
void ColoredDbg<qf_obj, key_obj>::copy(BitVector& src, BitVector& dst, uint64_t fromStartIdx, uint64_t toStartIdx) {
	uint64_t s = time(NULL);
	size_t i = 0;
    while (i < num_samples) {
      size_t bitCnt = std::min(num_samples-i, sizeof(size_t));
      size_t wrd = src.get_int(fromStartIdx+i, bitCnt);
      dst.set_int(toStartIdx+i, wrd, bitCnt);
      i+=bitCnt;
    }
	/*
	for (uint32_t i = 0; i < num_samples; i++, start_idx++)
		if (vector[i])
			bv_buffer.set(start_idx);
	*/
	eqt += time(NULL) - s;
}

template <class qf_obj, class key_obj>
void ColoredDbg<qf_obj, key_obj>::add_bitvector(BitVector& vector, uint64_t eq_id) {
	uint64_t s = time(NULL);
	uint64_t start_idx = (eq_id  % NUM_BV_BUFFER) * num_samples;
	size_t i = 0;
    while (i < num_samples) {
      size_t bitCnt = std::min(num_samples-i, sizeof(size_t));
      size_t wrd = vector.get_int(i, bitCnt);
      bv_buffer.set_int(start_idx+i, wrd, bitCnt);
      i+=bitCnt;
    }
	/*
	for (uint32_t i = 0; i < num_samples; i++, start_idx++)
		if (vector[i])
			bv_buffer.set(start_idx);
	*/
	eqt += time(NULL) - s;
}

template <class qf_obj, class key_obj>
void ColoredDbg<qf_obj, key_obj>::bv_buffer_serialize() {
	uint64_t s = time(NULL);
	accumEqClsCnt.push_back(get_num_eqclasses());
	BitVectorRRR final_com_bv(bv_buffer);
	rrrt += time(NULL) - s;
	std::string bv_file(prefix + std::to_string(num_serializations) + "_"
											EQCLASS_FILE);
	final_com_bv.serialize(bv_file);
	bv_buffer.reset();
	num_serializations++;
	eq_exportt += time(NULL) - s;
}

template <class qf_obj, class key_obj>
void ColoredDbg<qf_obj, key_obj>::serialize() {
	uint64_t s = time(NULL);
	// serialize the CQF
	dbg.serialize(prefix + CQF_FILE);

	cqf_exportt += time(NULL) - s;
	
	// serialize the bv buffer last time if needed
	if (get_num_eqclasses() % NUM_BV_BUFFER > 1)
		bv_buffer_serialize();
	
	//serialize the eq class id map
	std::ofstream opfile(prefix + SAMPLEID_FILE);
	for (auto sample : sampleid_map)
		opfile << sample.first << " " << sample.second << std::endl;
	opfile.close();

	//serialize the eq class accumulative counts per bucket
	std::ofstream aceqfile(prefix + ACCUM_EQCLS_CNT_FILE);
	for (auto cnt : accumEqClsCnt)
		aceqfile << cnt << "\n";
	aceqfile.close();
	
	// dump eq class abundance dist for further analysis.
	/*
	std::ofstream tmpfile(prefix + "eqclass_dist.lst");
	for (auto sample : eqclass_map)
		tmpfile << sample.second.first << " " << sample.second.second << std::endl;
	tmpfile.close();
	*/

	PRINT_CDBG("Timings:" << "\neq: " << eqt 
						  << "\ncqf: " << cqft
						  << "\nEq serializing: " << eq_exportt
						  << "\nRRR: " << rrrt
						  << "\ncqf serialization: " << cqf_exportt
						  << "\ncqf iteration and heap: " << cqf_heapt
						  << "\nmap+hash: " << map_hasht
						  << "\nreshuffle: " << reshufflet);
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
		uint64_t bucket_idx = start_idx / NUM_BV_BUFFER;
		uint64_t bucket_offset = (start_idx % NUM_BV_BUFFER) * num_samples;
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
	key_obj>::construct(qf_obj *incqfs, std::unordered_map<std::string,
											uint64_t>& cutoffs, cdbg_bv_map_t<__uint128_t,
											std::pair<uint64_t, uint64_t>>& map, uint64_t num_kmers)
{
	uint32_t nqf = num_samples;
	uint64_t counter = 0;
	uint64_t s = time(NULL);
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
			it_incqfs[i] = incqfs[i].obj->begin(it->second);
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
	cqf_heapt += time(NULL) - s;
	while (!minheap.empty()) {
		s = time(NULL);
		assert(minheap.size() == nqf);
		SampleObject<KeyObject> cur;
		BitVector eq_class(num_samples);
		// Get the smallest key from minheap and update the eqclass vector
		cur = minheap.top();
		minheap.pop();
		
		cqf_heapt += time(NULL) - s;
		
		s = time(NULL);
		eq_class.set(cur.id);
		eqt += time(NULL) - s;
		// Keep poping keys from minheap until you see a different key.
		// While poping keys build the eq class for cur.
		// Increment iterators for all CQFs whose keys are popped.
		while (!minheap.empty() && cur.obj.key == minheap.top().obj.key) {
			s = time(NULL);
			uint32_t id = minheap.top().id;
			minheap.pop();
			cqf_heapt += time(NULL) - s;
		
			s = time(NULL);
			eq_class.set(id);
			eqt += time(NULL) - s;
		
			s = time(NULL);
			++it_incqfs[id];
			if (it_incqfs[id].done())	// If the iterator is done then decrement nqf
				nqf--;
			else {	// Insert the current iterator head in minHeap
				KeyObject key = *it_incqfs[id];
				SampleObject<KeyObject> obj(key, incqfs[id].sample_id, id);
				minheap.push(obj);
			}
			cqf_heapt += time(NULL) - s;
		}
		s = time(NULL);
		// Move the iterator of the smallest key.
		++it_incqfs[cur.id];
		if (it_incqfs[cur.id].done())	// If the iterator is done then decrement nqf
			nqf--;
		else {	// Insert the current iterator head in minHeap
			KeyObject key = *it_incqfs[cur.id];
			SampleObject<KeyObject> obj(key, incqfs[cur.id].sample_id, cur.id);
			minheap.push(obj);
		}
		cqf_heapt += time(NULL) - s;
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
																				std::string& prefix, uint64_t nqf) :
	dbg(key_bits, seed), bv_buffer(NUM_BV_BUFFER * nqf),
	prefix(prefix), num_samples(nqf), num_serializations(0) {
		if (nqf < UINT64_MAX)
			is_sampling = true;
	}

template <class qf_obj, class key_obj>
ColoredDbg<qf_obj, key_obj>::ColoredDbg(std::string& cqf_file,
																				std::vector<std::string>&
																				eqclass_files, std::string&
																				sample_file)
: dbg(cqf_file, false), bv_buffer() {
	num_samples = 0;
	num_serializations = 0;

	std::map<int, std::string> sorted_files;
	for (std::string file : eqclass_files) {
		int id = std::stoi(first_part(last_part(file, '/'), '_'));
		sorted_files[id] = file;
	}

	for (auto file : sorted_files) {
		std::cout << "Reading eq class file: " << file.second << std::endl;
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
