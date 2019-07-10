/*
 * ============================================================================
 *
 *         Author:  Prashant Pandey (), ppandey@cs.stonybrook.edu
 *                  Mike Ferdman (), mferdman@cs.stonybrook.edu
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
#include "gqf_cpp.h"
#include "gqf/hashutil.h"
#include "common_types.h"
#include "mantisconfig.hpp"

#define MANTIS_DBG_IN_MEMORY (0x01)
#define MANTIS_DBG_ON_DISK (0x02)

typedef sdsl::bit_vector BitVector;
typedef sdsl::rrr_vector<63> BitVectorRRR;

struct hash128 {
	uint64_t operator()(const __uint128_t& val128) const
	{
		__uint128_t val = val128;
		// Using the same seed as we use in k-mer hashing.
		return MurmurHash64A((void*)&val, sizeof(__uint128_t),
												 2038074743);
	}
};

template <typename Key, typename Value>
using cdbg_bv_map_t = spp::sparse_hash_map<Key, Value, hash128>;

using default_cdbg_bv_map_t = cdbg_bv_map_t<__uint128_t,
			std::pair<uint64_t,uint64_t>>;

template <class qf_obj, class key_obj>
class ColoredDbg {
	public:
		ColoredDbg(std::string& cqf_file, std::vector<std::string>& eqclass_files,
							 std::string& sample_file, int flag);

		ColoredDbg(uint64_t qbits, uint64_t key_bits, enum qf_hashmode hashmode,
							 uint32_t seed, std::string& prefix, uint64_t nqf, int flag);

		void build_sampleid_map(qf_obj *incqfs);

		default_cdbg_bv_map_t&
			construct(qf_obj *incqfs, uint64_t num_kmers);

		void set_console(spdlog::logger* c) { console = c; }
		const CQF<key_obj> *get_cqf(void) const { return &dbg; }
		uint64_t get_num_bitvectors(void) const;
		uint64_t get_num_eqclasses(void) const { return eqclass_map.size(); }
		uint64_t get_num_samples(void) const { return num_samples; }
		std::string get_sample(uint32_t id) const;
		uint32_t seed(void) const { return dbg.seed(); }
		uint64_t range(void) const { return dbg.range(); }

		std::vector<uint64_t>
			find_samples(const mantis::QuerySet& kmers);

        std::unordered_map<uint64_t, std::vector<uint64_t>>
            find_samples(const std::unordered_map<mantis::KmerHash, uint64_t> &uniqueKmers);

		void serialize();
		void reinit(default_cdbg_bv_map_t& map);
		void set_flush_eqclass_dist(void) { flush_eqclass_dis = true; }







		// Mantis merge (Jamshed)

		// Required to instantitate the output CDBG.
		ColoredDbg(ColoredDbg<qf_obj, key_obj> &cdbg1, ColoredDbg<qf_obj, key_obj> &cdbg2,
					uint64_t qbits, std::string &prefix, int flag, std::string fileName = "");

		// Required to load the input CDBGs.
		ColoredDbg(std::string &cqfFile, std::string &sampleListFile,
					std::vector<std::string> &eqclassFiles, int flag);

		// Returns the vector of names of all the color(equivalence)-class bitvector files.
		std::vector<std::string> &get_eq_class_files();

		// Merges two Colored DBG objects dbg1 and dbg2 into this Colored DBG object.
		void construct(ColoredDbg<qf_obj, key_obj> &dbg1, ColoredDbg<qf_obj, key_obj> &dbg2);

		// Serialize the CQF and equivalence-class bitvectors.
		// TODO: Move the list concatenation phase inside the serialization to the constructor,
		// with introduction of a new field, just keep the disk-write part here; thus eliminating
		// the parameters.
		// TODO: After the earlier task, merge this overloaded method with the earlier one
		// through default argument.
		void serialize(ColoredDbg<qf_obj, key_obj> &dbg1, ColoredDbg<qf_obj, key_obj> &dbg2);






		
	private:
		// returns true if adding this k-mer increased the number of equivalence
		// classes
		// and false otherwise.
		bool add_kmer(const typename key_obj::kmer_t& hash, const BitVector&
									vector);
		void add_bitvector(const BitVector& vector, uint64_t eq_id);
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
		int dbg_alloc_flag;
		bool flush_eqclass_dis{false};
		std::time_t start_time_;
		spdlog::logger* console;







		// Mangtis merge: Jamshed

		
		// Equivalence / Color class bitvector file names for this DBG.
		std::vector<std::string> eqClsFiles;

		// Required for hashing pair objects for unordered_map.
		// TODO: Consult Professor on this hashing.
		struct PairHash
		{
			template<typename T1, typename T2>
			std::size_t operator() (const std::pair<T1, T2> &p) const
			{
				return std::hash<T1>()(p.first) ^ std::hash<T2>()(p.second);
			}
		};

		// 2-d grid of buckets, with bucket_(i, j) containing all equivalence class ID pairs that
		// read bitvectors from the i'th and the j'th block of the two to-be-merged DBGs, respectively.
		std::vector<std::vector<std::vector<std::pair<uint64_t, uint64_t>>>> bucket;

		// Hash map for equivalence class pairs in the format: < (id1, id2) --> (newID, abundance) >.
		std::unordered_map<std::pair<uint64_t, uint64_t>,
							std::pair<uint64_t, uint64_t>, PairHash> eqClsMap;

		// Gathers all the distinct equivalence ID pairs and their abundance from the two DBGs into
		// the hash map 'eqClsMap'. The new IDs are set to a default value of 0.
		void gather_distinct_eq_class_pairs(ColoredDbg<qf_obj, key_obj> &dbg1, ColoredDbg<qf_obj, key_obj> &dbg2);

		// Allocates space for the buckets.
		void init_bit_vec_block_buckets(ColoredDbg<qf_obj, key_obj> &dbg1, ColoredDbg<qf_obj, key_obj> &dbg2);

		// Insert the equivalence class pair (eqCls1, eqCls2) into the hash map 'eqClsMap'.
		void add_eq_class_pair(uint64_t eqCls1, uint64_t eqCls2);

		// Concatenates the BitVector objects bv1[id1] and bv2[id2] into resultVec;
		// where bv1 and bv2 are partial color class tables, with sample counts colCount1 and colCount2
		// respectively.
		void concat(const BitVectorRRR& bv1, const uint64_t colCount1, const uint64_t eqID1,
					const BitVectorRRR& bv2, const uint64_t colCount2, const uint64_t eqID2,
					BitVector &resultVec);

		// Concatenate required bit vectors from the DBGs into this merged DBG, to build the
		// merged color(equivalence)-class table.
		void build_eq_classes(ColoredDbg<qf_obj, key_obj> &dbg1, ColoredDbg<qf_obj, key_obj> &dbg2);

		// Merge the CQFs from the DBGs into the CQF of this DBG.
		void merge_CQFs(ColoredDbg<qf_obj, key_obj> &dbg1, ColoredDbg<qf_obj, key_obj> &dbg2);

		// Add a kmer with equivalence class ID 'eqID' into the CQF of this DBG.
		void add_kmer(uint64_t kmer, uint64_t eqID);

		// Returns the sample-id mapping.
		std::unordered_map<uint64_t, std::string> &get_sample_id_map();

		// Serialize the bitvector buffer to disk, when the current class count is 'eqClsCount'.
		// TODO: Merge this overloaded method with the earlier one through default argument.
		void bv_buffer_serialize(uint64_t eqClsCount);
};

template <class T>
class SampleObject {
	public:
		SampleObject() : obj(), sample_id(), id(0) {}
		SampleObject(T o, std::string& s = std::string(),
								 uint32_t id = 0) : obj(o), sample_id(s), id(id) {}
		SampleObject(const SampleObject& o) : obj(o.obj),
		sample_id(o.sample_id), id(o.id) {}

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
		total += eqclasses[i].size();

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
					 assert(it_local->second.first <= mantis::NUM_BV_BUFFER &&
									it_input.second.first <= mantis::NUM_BV_BUFFER);
					 uint64_t src_idx = ((it_local->second.first - 1) * num_samples);
					 uint64_t dest_idx = ((it_input.second.first - 1) * num_samples);
					 for (uint32_t i = 0; i < num_samples; i++, src_idx++, dest_idx++)
						 if (bv_buffer[src_idx])
							 new_bv_buffer[dest_idx] = 1;
				 }
			 }
			 bv_buffer = new_bv_buffer;
		 }

template <class qf_obj, class key_obj>
void ColoredDbg<qf_obj, key_obj>::reinit(cdbg_bv_map_t<__uint128_t,
																				 std::pair<uint64_t, uint64_t>>& map) {
	// dbg.reset();
	uint64_t qbits = log2(dbg.numslots());
	uint64_t keybits = dbg.keybits();
	enum qf_hashmode hashmode = dbg.hash_mode();
	uint64_t seed = dbg.seed();
	dbg.delete_file();
	CQF<key_obj>cqf(qbits, keybits, hashmode, seed, prefix + mantis::CQF_FILE);
	dbg = cqf;

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
bool ColoredDbg<qf_obj, key_obj>::add_kmer(const typename key_obj::kmer_t&
																					 key, const BitVector& vector) {
	// A kmer (hash) is seen only once during the merge process.
	// So we insert every kmer in the dbg
	uint64_t eq_id;
	__uint128_t vec_hash = MurmurHash128A((void*)vector.data(),
																				vector.capacity()/8, 2038074743,
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

	// check: the k-mer should not already be present.
	uint64_t count = dbg.query(KeyObject(key,0,eq_id), QF_NO_LOCK |
													   QF_KEY_IS_HASH);
	if (count > 0) {
		console->error("K-mer was already present. kmer: {} eqid: {}", key, count);
		exit(1);
	}

	// we use the count to store the eqclass ids
	int ret = dbg.insert(KeyObject(key,0,eq_id), QF_NO_LOCK | QF_KEY_IS_HASH);
	if (ret == QF_NO_SPACE) {
		// This means that auto_resize failed. 
		console->error("The CQF is full and auto resize failed. Please rerun build with a bigger size.");
		exit(1);
	}

	return added_eq_class;
}

template <class qf_obj, class key_obj>
void ColoredDbg<qf_obj, key_obj>::add_bitvector(const BitVector& vector,
																								uint64_t eq_id) {
	uint64_t start_idx = (eq_id  % mantis::NUM_BV_BUFFER) * num_samples;
	for (uint32_t i = 0; i < num_samples/64*64; i+=64)
		bv_buffer.set_int(start_idx+i, vector.get_int(i, 64), 64);
	if (num_samples%64)
		bv_buffer.set_int(start_idx+num_samples/64*64,
											vector.get_int(num_samples/64*64, num_samples%64),
											num_samples%64);
}

template <class qf_obj, class key_obj>
void ColoredDbg<qf_obj, key_obj>::bv_buffer_serialize() {
	BitVector bv_temp(bv_buffer);
	if (get_num_eqclasses() % mantis::NUM_BV_BUFFER > 0) {
		bv_temp.resize((get_num_eqclasses() % mantis::NUM_BV_BUFFER) *
									 num_samples);
	}

	BitVectorRRR final_com_bv(bv_temp);
	std::string bv_file(prefix + std::to_string(num_serializations) + "_" +
											mantis::EQCLASS_FILE);
	sdsl::store_to_file(final_com_bv, bv_file);
	bv_buffer = BitVector(bv_buffer.bit_size());
	num_serializations++;
}

template <class qf_obj, class key_obj>
void ColoredDbg<qf_obj, key_obj>::serialize() {
	// serialize the CQF
	if (dbg_alloc_flag == MANTIS_DBG_IN_MEMORY)
		dbg.serialize(prefix + mantis::CQF_FILE);
	else
		dbg.close();

	// serialize the bv buffer last time if needed
	if (get_num_eqclasses() % mantis::NUM_BV_BUFFER > 0)
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
			tmpfile << sample.second.first << " " << sample.second.second <<
				std::endl;
		tmpfile.close();
	}
}

template <class qf_obj, class key_obj>
std::vector<uint64_t>
ColoredDbg<qf_obj,key_obj>::find_samples(const mantis::QuerySet& kmers) {
	// Find a list of eq classes and the number of kmers that belong those eq
	// classes.
	std::unordered_map<uint64_t, uint64_t> query_eqclass_map;
	for (auto k : kmers) {
		key_obj key(k, 0, 0);
		uint64_t eqclass = dbg.query(key, 0);
		if (eqclass)
			query_eqclass_map[eqclass] += 1;
	}

	std::vector<uint64_t> sample_map(num_samples, 0);
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
std::unordered_map<uint64_t, std::vector<uint64_t>>
ColoredDbg<qf_obj,key_obj>::find_samples(const std::unordered_map<mantis::KmerHash, uint64_t> &uniqueKmers) {
	// Find a list of eq classes and the number of kmers that belong those eq
	// classes.
	std::unordered_map<uint64_t, std::vector<uint64_t>> query_eqclass_map;
	for (auto kv : uniqueKmers) {
		key_obj key(kv.first, 0, 0);
		uint64_t eqclass = dbg.query(key, 0);
		if (eqclass) {
		    kv.second = eqclass;
            query_eqclass_map[eqclass] = std::vector<uint64_t>();
        }
	}

	std::vector<uint64_t> sample_map(num_samples, 0);
	for (auto it = query_eqclass_map.begin(); it != query_eqclass_map.end();
		 ++it) {
		auto eqclass_id = it->first;
		auto &vec = it->second;
		// counter starts from 1.
		uint64_t start_idx = (eqclass_id - 1);
		uint64_t bucket_idx = start_idx / mantis::NUM_BV_BUFFER;
		uint64_t bucket_offset = (start_idx % mantis::NUM_BV_BUFFER) * num_samples;
		for (uint32_t w = 0; w <= num_samples / 64; w++) {
			uint64_t len = std::min((uint64_t)64, num_samples - w * 64);
			uint64_t wrd = eqclasses[bucket_idx].get_int(bucket_offset, len);
			for (uint32_t i = 0, sCntr = w * 64; i < len; i++, sCntr++)
				if ((wrd >> i) & 0x01)
					vec.push_back(sCntr);
			bucket_offset += len;
		}
	}
	return query_eqclass_map;
}

template <class qf_obj, class key_obj>
cdbg_bv_map_t<__uint128_t, std::pair<uint64_t, uint64_t>>& ColoredDbg<qf_obj,
	key_obj>::construct(qf_obj *incqfs, uint64_t num_kmers)
{
	uint64_t counter = 0;
	bool is_sampling = (num_kmers < std::numeric_limits<uint64_t>::max());

	struct Iterator {
		QFi qfi;
		typename key_obj::kmer_t kmer;
		uint32_t id;
		bool do_madvice{false};
		Iterator(uint32_t id, const QF* cqf, bool flag): id(id), do_madvice(flag)
		{
			if (qf_iterator_from_position(cqf, &qfi, 0) != QFI_INVALID) {
				get_key();
        if (do_madvice)
          qfi_initial_madvise(&qfi);
      }
		}
		bool next() {
			if (do_madvice) {
				if (qfi_next_madvise(&qfi) == QFI_INVALID) return false;
			} else {
				if (qfi_next(&qfi) == QFI_INVALID) return false;
			}
			get_key();
			return true;
		}
		bool end() const {
			return qfi_end(&qfi);
		}
		bool operator>(const Iterator& rhs) const {
			return key() > rhs.key();
		}
		const typename key_obj::kmer_t& key() const { return kmer; }
		private:
		void get_key() {
			uint64_t value, count;
			qfi_get_hash(&qfi, &kmer, &value, &count);
		}
	};

  typename CQF<key_obj>::Iterator walk_behind_iterator;
  
	struct Minheap_PQ {
		void push(const Iterator& obj) {
			c.emplace_back(obj);
			std::push_heap(c.begin(), c.end(), std::greater<Iterator>());
		}
		void pop() {
			std::pop_heap(c.begin(), c.end(), std::greater<Iterator>());
			c.pop_back();
		}
		void replace_top(const Iterator& obj) {
			c.emplace_back(obj);
			pop();
		}
		Iterator& top() { return c.front(); }
		bool empty() const { return c.empty(); }
		private:
		std::vector<Iterator> c;
	};
	Minheap_PQ minheap;

	for (uint32_t i = 0; i < num_samples; i++) {
		Iterator qfi(i, incqfs[i].obj->get_cqf(), true);
		if (qfi.end()) continue;
		minheap.push(qfi);
	}

	while (!minheap.empty()) {
		BitVector eq_class(num_samples);
		KeyObject::kmer_t last_key;
		do {
			Iterator& cur = minheap.top();
			last_key = cur.key();
			eq_class[cur.id] = 1;
			if (cur.next())
				minheap.replace_top(cur);
			else
				minheap.pop();
		} while(!minheap.empty() && last_key == minheap.top().key());
		bool added_eq_class = add_kmer(last_key, eq_class);
		++counter;

    if (counter == 4096) {
      walk_behind_iterator = dbg.begin();
    } else if (counter > 4096) {
      ++walk_behind_iterator;
    }
    
		// Progress tracker
		static uint64_t last_size = 0;
		if (dbg.dist_elts() % 10000000 == 0 &&
				dbg.dist_elts() != last_size) {
			last_size = dbg.dist_elts();
			console->info("Kmers merged: {}  Num eq classes: {}  Total time: {}",
										dbg.dist_elts(), get_num_eqclasses(), time(nullptr) -
										start_time_);
		}

		// Check if the bit vector buffer is full and needs to be serialized.
		if (added_eq_class and (get_num_eqclasses() % mantis::NUM_BV_BUFFER == 0))
		{
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

		//while(!minheap.empty() && minheap.top().end()) minheap.pop();
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
																				enum qf_hashmode hashmode,
																				uint32_t seed, std::string& prefix,
																				uint64_t nqf, int flag) :
	bv_buffer(mantis::NUM_BV_BUFFER * nqf), prefix(prefix), num_samples(nqf),
	num_serializations(0), start_time_(std::time(nullptr)) {
		if (flag == MANTIS_DBG_IN_MEMORY) {
			CQF<key_obj> cqf(qbits, key_bits, hashmode, seed);
			dbg = cqf;
			dbg_alloc_flag = MANTIS_DBG_IN_MEMORY;
		}
		else if (flag == MANTIS_DBG_ON_DISK) {
			CQF<key_obj> cqf(qbits, key_bits, hashmode, seed, prefix +
                       mantis::CQF_FILE);
			dbg = cqf;
			dbg_alloc_flag = MANTIS_DBG_ON_DISK;
		} else {
			ERROR("Wrong Mantis alloc mode.");
			exit(EXIT_FAILURE);
		}
		dbg.set_auto_resize();
	}

template <class qf_obj, class key_obj>
ColoredDbg<qf_obj, key_obj>::ColoredDbg(std::string& cqf_file,
																				std::vector<std::string>&
																				eqclass_files, std::string&
																				sample_file, int flag) : bv_buffer(),
	start_time_(std::time(nullptr)) {
		num_samples = 0;
		num_serializations = 0;

		if (flag == MANTIS_DBG_IN_MEMORY) {
			CQF<key_obj>cqf(cqf_file, CQF_FREAD);
			dbg = cqf;
			dbg_alloc_flag = MANTIS_DBG_IN_MEMORY;
		} else if (flag == MANTIS_DBG_ON_DISK) {
			CQF<key_obj>cqf(cqf_file, CQF_MMAP);
			dbg = cqf;
			dbg_alloc_flag = MANTIS_DBG_ON_DISK;
		} else {
			ERROR("Wrong Mantis alloc mode.");
			exit(EXIT_FAILURE);
		}

		std::map<int, std::string> sorted_files;
		for (std::string file : eqclass_files) {
			int id = std::stoi(first_part(last_part(file, '/'), '_'));
			sorted_files[id] = file;
		}

		eqclasses.reserve(sorted_files.size());
		BitVectorRRR bv;
		for (auto file : sorted_files) {
			sdsl::load_from_file(bv, file.second);
			eqclasses.push_back(bv);
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







// Mantis merge: Jamshed

template <class qf_obj, class key_obj>
ColoredDbg<qf_obj, key_obj>::ColoredDbg(ColoredDbg<qf_obj, key_obj> &cdbg1, ColoredDbg<qf_obj, key_obj> &cdbg2,
										uint64_t qbits, std::string &prefix, int flag, std::string fileName):
	bv_buffer(mantis::NUM_BV_BUFFER * (cdbg1.get_num_samples() + cdbg2.get_num_samples())),
	prefix(prefix), num_samples(cdbg1.get_num_samples() + cdbg2.get_num_samples()),
	num_serializations(0), start_time_(std::time(nullptr))
	{
		if(flag == MANTIS_DBG_IN_MEMORY)
		{
			CQF<key_obj> cqf(qbits, cdbg1.get_cqf() -> keybits(), cdbg1.get_cqf() -> hash_mode(),
							cdbg1.get_cqf() -> seed());
			dbg = cqf;
			dbg_alloc_flag = MANTIS_DBG_IN_MEMORY;
		}
		else if(flag == MANTIS_DBG_ON_DISK)
		{
			CQF<key_obj> cqf(qbits, cdbg1.get_cqf() -> keybits(), cdbg1.get_cqf() -> hash_mode(),
							cdbg1.get_cqf() -> seed(), prefix + fileName);
			dbg = cqf;
			dbg_alloc_flag = MANTIS_DBG_ON_DISK;
		}
		else
		{
			ERROR("Wrong Mantis alloc mode.");
			exit(EXIT_FAILURE);
		}

		dbg.set_auto_resize();

		// printf("\n\nMSG: CDBG (output) constructed; initial size = %llu.\n\n",
		// 		(unsigned long long)dbg.dist_elts());
	}





template <class qf_obj, class key_obj>
ColoredDbg<qf_obj, key_obj> ::
	ColoredDbg(std::string &cqfFile, std::string &sampleListFile, std::vector<std::string> &eqclassFiles,
				int flag):
	bv_buffer(), start_time_(std::time(nullptr))
	{
		num_samples = 0;
		num_serializations = 0;


		// Load the CQF
		if(flag == MANTIS_DBG_IN_MEMORY)
		{
			CQF<key_obj> cqf(cqfFile, CQF_FREAD);
			dbg = cqf;
			dbg_alloc_flag = MANTIS_DBG_IN_MEMORY;
		}
		else if(flag == MANTIS_DBG_ON_DISK)
		{
			CQF<key_obj> cqf(cqfFile, CQF_MMAP);
			dbg = cqf;
			dbg_alloc_flag = MANTIS_DBG_ON_DISK;
		}
		else
		{
			ERROR("Wrong Mantis alloc mode.");
			exit(EXIT_FAILURE);
		}


		//printf("\n\nLoaded CQF size = %llu\n\n", (unsigned long long)dbg.dist_elts());


		// Load the sample / experiment names.
		std::ifstream sampleList(sampleListFile.c_str());
		std::string sampleName;
		uint32_t sampleID;

		while (sampleList >> sampleID >> sampleName)
		{
			sampleid_map[sampleID] = sampleName;
			num_samples++;
		}

		sampleList.close();


		// Load the equivalence class bitvector file names only;
		// not loading the bitvectors (color classes) all at once, like the other identical constructor.
		std::map<int, std::string> sortedFiles;
		for (std::string file : eqclassFiles)
		{
			int fileID = std::stoi(first_part(last_part(file, '/'), '_'));
			sortedFiles[fileID] = file;

			printf("\nLoaded color-class file %s\n", file.c_str());
		}
		

		eqClsFiles.reserve(eqclassFiles.size());
		for(auto idFilePair : sortedFiles)
			eqClsFiles.push_back(idFilePair.second);

		num_serializations = eqClsFiles.size();


		/*
		eqclasses.reserve(sorted_files.size());
		BitVectorRRR bv;
		for (auto file : sorted_files) {
			sdsl::load_from_file(bv, file.second);
			eqclasses.push_back(bv);
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
		*/
}





template <typename qf_obj, typename key_obj>
inline std::vector<std::string> &ColoredDbg<qf_obj, key_obj> ::
	get_eq_class_files()
{
	return eqClsFiles;
}





template <typename qf_obj, typename key_obj>
inline void ColoredDbg<qf_obj, key_obj> ::
	init_bit_vec_block_buckets(ColoredDbg<qf_obj, key_obj> &dbg1, ColoredDbg<qf_obj, key_obj> &dbg2)
{
	const uint64_t fileCount1 = (dbg1.get_num_eqclasses() / mantis::NUM_BV_BUFFER) + 1,
					fileCount2 = (dbg2.get_num_eqclasses() / mantis::NUM_BV_BUFFER) + 1;
	
	bucket.resize(fileCount1);
	for(uint64_t i = 0; i < fileCount1; ++i)
		bucket[i].resize(fileCount2);
}





template <typename qf_obj, typename key_obj>
inline void ColoredDbg<qf_obj, key_obj> ::
	add_eq_class_pair(uint64_t eqCls1, uint64_t eqCls2)
{
	auto eqPair = std::make_pair(eqCls1, eqCls2);
	auto it = eqClsMap.find(eqPair);

	const static auto newEntry = std::make_pair((uint64_t)0, (uint64_t)1);

	if(it == eqClsMap.end())
	{
		eqClsMap[eqPair] = newEntry;

		bucket[eqCls1 / mantis::NUM_BV_BUFFER][eqCls2 / mantis::NUM_BV_BUFFER].push_back(eqPair);

		printf("New equivalence pair: (%d, %d)\n", (int)eqPair.first, (int)eqPair.second);
	}
	else
		it -> second.second++;
}





template <typename qf_obj, typename key_obj>
void ColoredDbg<qf_obj, key_obj> ::
	gather_distinct_eq_class_pairs(ColoredDbg<qf_obj, key_obj> &dbg1, ColoredDbg<qf_obj, key_obj> &dbg2)
{
	puts("\nMSG: At distinct equivalence-class pair gathering phase.\n");

	const CQF<key_obj> *cqf1 = dbg1.get_cqf(), *cqf2 = dbg2.get_cqf();
	CQF<KeyObject> :: Iterator it1 = cqf1 -> begin(), it2 = cqf2 -> begin();

	uint64_t kmerCount = 0, equalKmerCount = 0;	// for debugging purposes

	
	init_bit_vec_block_buckets(dbg1, dbg2);


	if(it1.done() || it2.done())
		puts("\n\nERR MSG: One or more iterator already at end position before starting walk.\n\n");

	printf("\nDBG1 size = %d, DBG2 size = %d\n", (int)cqf1 -> dist_elts(), (int)cqf2 -> dist_elts());


	while(!it1.done() || !it2.done())
	{
		uint64_t kmer1, kmer2, eqClass1, eqClass2;
		key_obj cqfEntry1, cqfEntry2;

		if(!it1.done())
		{
			cqfEntry1 = it1.get_cur_hash();
			kmer1 = cqfEntry1.key;
		}
		
		if(!it2.done())
		{
			cqfEntry2 = it2.get_cur_hash();
			kmer2 = cqfEntry2.key;
		}


		// eqClassX = 0 implies absence in CDBG X.


		if(it1.done())
			eqClass1 = 0, eqClass2 = cqfEntry2.count, ++it2;
		else if(it2.done())
			eqClass1 = cqfEntry1.count, eqClass2 = 0, ++it1;
		else if(kmer1 < kmer2)
			eqClass1 = cqfEntry1.count, eqClass2 = 0, ++it1;
		else if(kmer2 < kmer1)
			eqClass1 = 0, eqClass2 = cqfEntry2.count, ++it2;
		else
			eqClass1 = cqfEntry1.count, eqClass2 = cqfEntry2.count, ++it1, ++it2,
			equalKmerCount++;


		add_eq_class_pair(eqClass1, eqClass2);
		kmerCount++;
	}


	printf("\n\nMSG: Distinct kmers found %d, equal kmers found %d (at distinct pairs gathering phase).\n\n",
			(int)kmerCount, (int32_t)equalKmerCount);
}





template <typename qf_obj, typename key_obj>
void ColoredDbg<qf_obj, key_obj> ::
	concat(const BitVectorRRR& bv1, const uint64_t colCount1, const uint64_t eqID1,
			const BitVectorRRR& bv2, const uint64_t colCount2, const uint64_t eqID2,
			BitVector &resultVec)
{
	const uint64_t wordLen = 64;


	if(eqID1)
	{
		uint64_t offset1 = ((eqID1 - 1) % mantis::NUM_BV_BUFFER) * colCount1;

		for(uint32_t wordCount = 0; wordCount <= colCount1 / wordLen; ++wordCount)
		{
			uint64_t readLen = std::min(wordLen, colCount1 - wordCount * wordLen);
			uint64_t word = bv1.get_int(offset1, readLen);

			// Optimize here; preferrably eliminate the loop with one statement (some sort of set_int() ?).
			for(uint32_t bitIdx = 0, sampleCounter = wordCount * wordLen; bitIdx < readLen; ++bitIdx, ++sampleCounter)
				if((word >> bitIdx) & 0x01)
					resultVec[sampleCounter] = 1;				

			offset1 += readLen;
		}
	}


	if(eqID2)
	{
		u_int64_t offset2 = ((eqID2 - 1) % mantis::NUM_BV_BUFFER) * colCount2;

		for(uint32_t wordCount = 0; wordCount <= colCount2 / wordLen; ++wordCount)
		{
			uint64_t readLen = std::min(wordLen, colCount2 - wordCount * wordLen);
			uint64_t word = bv2.get_int(offset2, readLen);

			// Optimize here; preferrably eliminate the loop with one statement (some sort of set_int() ?).
			for(uint32_t bitIdx = 0, sampleCounter = wordCount * wordLen; bitIdx < readLen; ++bitIdx, ++sampleCounter)
				if((word >> bitIdx) & 0x01)
					resultVec[colCount1 + sampleCounter] = 1;

			offset2 += readLen;
		}
	}

	
	//std::vector<uint64_t> sample_map(num_samples, 0);
	//for (auto it = query_eqclass_map.begin(); it != query_eqclass_map.end(); ++it) {
		// auto eqclass_id = it->first;
		// auto count = it->second;
		// counter starts from 1.
		// uint64_t start_idx = (eqclass_id - 1);
		// uint64_t bucket_idx = start_idx / mantis::NUM_BV_BUFFER;
		// uint64_t bucket_offset = (start_idx % mantis::NUM_BV_BUFFER) * num_samples;
		// for (uint32_t w = 0; w <= num_samples / 64; w++) {
		// 	uint64_t len = std::min((uint64_t)64, num_samples - w * 64);
		// 	uint64_t wrd = eqclasses[bucket_idx].get_int(bucket_offset, len);
		// 	for (uint32_t i = 0, sCntr = w * 64; i < len; i++, sCntr++)
		// 		if ((wrd >> i) & 0x01)
		// 			sample_map[sCntr] += count;
		// 	bucket_offset += len;
	//	}
}





template <class qf_obj, class key_obj>
void ColoredDbg<qf_obj, key_obj>::bv_buffer_serialize(uint64_t eqClsCount)
{
	BitVector bv_temp(bv_buffer);
	if (eqClsCount % mantis::NUM_BV_BUFFER > 0)
		bv_temp.resize((eqClsCount % mantis::NUM_BV_BUFFER) * num_samples);

	const char OUTPUT_EQCLASS_FILE[] = "merged_eqclass_rrr.cls";

	BitVectorRRR final_com_bv(bv_temp);
	//std::string bv_file(prefix + std::to_string(num_serializations) + "_" + mantis::EQCLASS_FILE);
	std::string bv_file(prefix + std::to_string(num_serializations) + "_" + OUTPUT_EQCLASS_FILE);

	sdsl::store_to_file(final_com_bv, bv_file);
	bv_buffer = BitVector(bv_buffer.bit_size());

	num_serializations++;
}





template <typename qf_obj, typename key_obj>
void ColoredDbg<qf_obj, key_obj> ::
	build_eq_classes(ColoredDbg<qf_obj, key_obj> &dbg1, ColoredDbg<qf_obj, key_obj> &dbg2)
{
	puts("\nMSG: At color-class building (bitvectors concatenation) phase.\n");

	const uint64_t fileCount1 = (dbg1.get_num_eqclasses() / mantis::NUM_BV_BUFFER) + 1,
					fileCount2 = (dbg2.get_num_eqclasses() / mantis::NUM_BV_BUFFER) + 1;
	uint64_t serialID = 0;


	// No need for the following statement. eqclass_files should come with the constructor of the dbgs
	// eqclass_files = mantis::fs::GetFilesExt(prefix.c_str(), mantis::EQCLASS_FILE);

	// No need of sorting; done at new constructor.
	// required sort of the file names (?)


	for(uint64_t i = 0; i < fileCount1; ++i)
	{
		// Required: data_1 = read_i'th Bit Vector block for dbg1
		BitVectorRRR bitVec1;
		sdsl::load_from_file(bitVec1, dbg1.get_eq_class_files()[i]);

		for(uint64_t j = 0; j < fileCount2; ++j)
		{
			printf("\nAt bucket (%d, %d). Size = %d\n", (int)i, (int)j, (int)bucket[i][j].size());
			

			// Required: data_2 = read_j'th Bit Vector block for dbg2
			BitVectorRRR bitVec2;
			sdsl::load_from_file(bitVec2, dbg2.get_eq_class_files()[j]);

			for(auto eqClsPair : bucket[i][j])
			{
				uint64_t eq1 = eqClsPair.first, eq2 = eqClsPair.second;
				eqClsMap[eqClsPair].first = ++serialID; // serialID++ doesn't work. Introduces bug in CQF;
														// Probable cause: CQF doesn't support 0-count ?

				// Required: bit_vector = concat(data_1.query(eq1), data_2.query(eq2))
				// Required: dbg.bitVectorBuffer[serialID] = bit_vector
				

				BitVector mergedEqCls(num_samples);

				concat(bitVec1, dbg1.get_num_samples(), eq1, bitVec2, dbg2.get_num_samples(), eq2, mergedEqCls);
				add_bitvector(mergedEqCls, serialID);

				
				// Serialization and disk-write if required
				if(serialID % mantis::NUM_BV_BUFFER == 0)
					bv_buffer_serialize(serialID);
			}
		}
	}


	// Serialize the bv buffer last time if needed
	if (serialID % mantis::NUM_BV_BUFFER > 0)
		bv_buffer_serialize(serialID);

	
	// puts("\nMSG: Done building.\n");
}





template <typename qf_obj, typename key_obj>
inline void ColoredDbg<qf_obj, key_obj> ::
	add_kmer(uint64_t kmer, uint64_t eqID)
{
	if (dbg.insert(KeyObject(kmer, 0, eqID), QF_NO_LOCK | QF_KEY_IS_HASH) == QF_NO_SPACE)
	{
		// This means that auto_resize failed.
		console -> error("The CQF is full and auto resize failed. Please rerun build with a bigger size.");
		exit(1);
	}
}





template <typename qf_obj, typename key_obj>
void ColoredDbg<qf_obj, key_obj> ::
	merge_CQFs(ColoredDbg<qf_obj, key_obj> &dbg1, ColoredDbg<qf_obj, key_obj> &dbg2)
{
	puts("\nMSG: At CQFs merging phase.\n");

	const CQF<key_obj> *cqf1 = dbg1.get_cqf(), *cqf2 = dbg2.get_cqf();
	CQF<KeyObject> :: Iterator it1 = cqf1 -> begin(), it2 = cqf2 -> begin();

	uint64_t kmerCount = 0, equalKmerCount = 0; // for debugging purpose(s)

	
	if(it1.done() || it2.done())
		puts("\n\nMSG: One or more iterator already at end position before starting walk.\n\n");


	while(!it1.done() || !it2.done())
	{
		uint64_t kmer1, kmer2, kmer, eqClass1, eqClass2;
		key_obj cqfEntry1, cqfEntry2;

		if(!it1.done())
		{
			cqfEntry1 = it1.get_cur_hash();
			kmer1 = cqfEntry1.key;
		}
		
		if(!it2.done())
		{
			cqfEntry2 = it2.get_cur_hash();
			kmer2 = cqfEntry2.key;
		}


		if(it1.done())
		{
			kmer = kmer2;
			eqClass1 = 0, eqClass2 = cqfEntry2.count;

			++it2;
		}
		else if(it2.done())
		{
			kmer = kmer1;
			eqClass1 = cqfEntry1.count, eqClass2 = 0;

			++it1;
		}
		else if(kmer1 < kmer2)
		{
			kmer = kmer1;
			eqClass1 = cqfEntry1.count, eqClass2 = 0;

			++it1;
		}
		else if(kmer2 < kmer1)
		{
			kmer = kmer2;
			eqClass1 = 0, eqClass2 = cqfEntry2.count;

			++it2;
		}
		else
		{
			kmer = kmer1;
			eqClass1 = cqfEntry1.count, eqClass2 = cqfEntry2.count;

			++it1, ++it2;

			equalKmerCount++; // for debugging purpose(s)
		}


		add_kmer(kmer, eqClsMap[std :: make_pair(eqClass1, eqClass2)].first);
		kmerCount++;
	}


	printf("\nMSG: In the merged CQF: #kmers = %llu; #equal_kmers = %llu\n\n",
			(unsigned long long)kmerCount, (unsigned long long)equalKmerCount);
	printf("\nMSG: Merged CQF size: %llu\n", (unsigned long long)dbg.dist_elts());
}





template <typename qf_obj, typename key_obj>
void ColoredDbg<qf_obj, key_obj> :: 
	construct(ColoredDbg<qf_obj, key_obj> &dbg1, ColoredDbg<qf_obj, key_obj> &dbg2)
{
	puts("\n\nMSG: Merge starting.\n\n");

	gather_distinct_eq_class_pairs(dbg1, dbg2);

	build_eq_classes(dbg1, dbg2);

	merge_CQFs(dbg1, dbg2);

	printf("\n\nSize of merged CQF = %llu\n\n", (unsigned long long)dbg.dist_elts());

	puts("\n\nMerge ending.\n\n");
}





template <typename qf_obj, typename key_obj>
std::unordered_map<uint64_t, std::string> &ColoredDbg<qf_obj, key_obj>::get_sample_id_map()
{
	return sampleid_map;
}





template <typename qf_obj, typename key_obj>
void ColoredDbg<qf_obj, key_obj>::
	serialize(ColoredDbg<qf_obj, key_obj> &dbg1, ColoredDbg<qf_obj, key_obj> &dbg2)
{
	// Serialize the CQF
	const char OUTPUT_CQF_FILE[] = "merged_dbg_cqf.ser";

	if (dbg_alloc_flag == MANTIS_DBG_IN_MEMORY)
		//dbg.serialize(prefix + mantis::CQF_FILE);
		dbg.serialize(prefix + OUTPUT_CQF_FILE);
	else
		dbg.close();


	// Serialize the bv buffer last time if needed.
	// Done at the color(equivalence)-class building phase.


	// TODO: Move the list concatenation phase to constructor, with introduction of a new field;
	// just keep the disk-write part here.

	// Serialize the eq class id map
	char OUTPUT_SAMPLE_ID_LIST[] = "merged_sampleid.lst";
	//std::ofstream outputFile(prefix + mantis::SAMPLEID_FILE);
	std::ofstream outputFile(prefix + OUTPUT_SAMPLE_ID_LIST);

	for(auto sampleID: dbg1.get_sample_id_map())
		outputFile << sampleID.first << " " << sampleID.second << std::endl;

	for(auto sampleID: dbg2.get_sample_id_map())
		outputFile << dbg1.get_num_samples() + sampleID.first << " " << sampleID.second << std::endl;

	outputFile.close();


	// Dump the abundance distribution of the equivalence / color classes.
	if (flush_eqclass_dis)
	{
		const char OUTPUT_ABUNDANCE_DIST_FILE[] = "merged_eqclass_dist.lst";
		std::ofstream outputFile(prefix + OUTPUT_ABUNDANCE_DIST_FILE);

		for (auto idFreq : eqClsMap)
			outputFile << idFreq.second.first << " " << idFreq.second.second << std::endl;
		
		outputFile.close();
	}
}

#endif
