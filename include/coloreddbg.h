/*
 * ============================================================================
 *
 *         Author:  Prashant Pandey (), ppandey@cs.stonybrook.edu
 *                  Mike Ferdman (), mferdman@cs.stonybrook.edu
 * 					Jamshed Khan (), jamshed@cs.umd.edu
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

#include "BooPHF.h"

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
					std::string &prefix, int flag);

		// Required to load the input CDBGs.
		ColoredDbg(std::string &cqfFile, std::string &sampleListFile,
					std::vector<std::string> &eqclassFiles, int flag);

		// Returns the vector of names of all the color(equivalence)-class bitvector files.
		inline std::vector<std::string> &get_eq_class_files() { return eqClsFiles; }

		// Returns the number of color(equivalence)-class bitvector files.
		inline uint64_t get_eq_class_file_count() { return eqClsFiles.size(); }

		// Returns the number of equivalence / color classes in the merged CdBG.
		uint64_t get_eqclass_count()	{ return eqClsMap.size(); }

		// Merges two Colored DBG objects dbg1 and dbg2 into this Colored DBG object.
		void construct(ColoredDbg<qf_obj, key_obj> &dbg1, ColoredDbg<qf_obj, key_obj> &dbg2);

		// Serialize the CQF and equivalence-class bitvectors.
		// TODO: Move the list concatenation phase inside the serialization to the constructor,
		// with introduction of a new field, just keep the disk-write part here; thus eliminating
		// the parameters.
		// TODO: After the earlier task, merge this overloaded method with the earlier one
		// through default argument.
		void serialize(ColoredDbg<qf_obj, key_obj> &dbg1, ColoredDbg<qf_obj, key_obj> &dbg2);

		// Returns the collection of BitvectorRRR's (compressed color-classes) of this CdBG.
		std::vector<BitVectorRRR> get_eqclasses() { return eqclasses; }





		// Merge approach 2 (Jamshed)

		// Merges two Colored dBG objects 'cdbg1' and 'cdbg2' into this Colored dBG;
		// returns the number of distinct color-classes at the merged CdBG.
		uint64_t merge(ColoredDbg<qf_obj, key_obj> &cdbg1, ColoredDbg<qf_obj, key_obj> &cdbg2);

		// Sets the number of processor-threads to be used at the intermediate steps of 
		// unique id-pairs filtering and MPH building.
		inline void set_thread_count(uint threadNum) { threadCount = threadNum; }

		// Sets the maximum memory usage limit for the intermediate step of unique id-pairs filtering.
		inline void set_max_memory_for_sort(uint maxMemory) { maxMemoryForSort = maxMemory; }

		// Serializes the output CQF, sample-id mapping, and abundance distribution (optional).
		void serialize_cqf_and_abundance_dist(ColoredDbg<qf_obj, key_obj> &cdbg1, ColoredDbg<qf_obj, key_obj> &cbg2);





		// Merge approach 3 (Jamshed)

		// Merges two Colored dBG objects 'cdbg1' and 'cdbg2' into this Colored dBG;
		// returns the number of distinct color-classes at the merged CdBG.
		uint64_t merge_2(ColoredDbg<qf_obj, key_obj> &cdbg1, ColoredDbg<qf_obj, key_obj> &cdbg2);




		
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

		const static uint64_t PROGRESS_STEP = 10000000;
		const static uint64_t ITERATOR_WINDOW_SIZE = 4096;

		
		// Equivalence / Color class bitvector file names for this DBG.
		std::vector<std::string> eqClsFiles;

		// Required for hashing pair objects for unordered_map.
		// TODO: Consult Professor on this hashing.
		// Done. Resorted to boost::hash_combine from plain XOR hashing.
		// For more explanation, consult https://stackoverflow.com/questions/35985960/c-why-is-boosthash-combine-the-best-way-to-combine-hash-values
		struct PairHash
		{
			template<typename T1, typename T2>
			std::size_t operator() (const std::pair<T1, T2> &p) const
			{
				size_t seed = std::hash<T1>()(p.first) + 0x9e3779b9;
				seed ^= std::hash<T1>()(p.second) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
				return seed;
			}
		};

		// 2-d grid of buckets, with bucket_(i, j) containing all equivalence class ID pairs that
		// read bitvectors from the i'th and the j'th block of the two to-be-merged DBGs, respectively.
		std::vector<std::vector<std::vector<std::pair<uint64_t, uint64_t>>>> bucket;

		// Hash map for equivalence class pairs in the format: < (id1, id2) --> (newID, abundance) >.
		std::unordered_map<std::pair<uint64_t, uint64_t>,
							std::pair<uint64_t, uint64_t>, PairHash> eqClsMap;

		// Advances the CQF iterator 'it', with keeping track of the 'step' count; and fetches the next
		// CQF-entry into 'cqfEntry' if the iterator 'it' is advanced into a non-end position. Also, advances
		// or initializes the iterator 'walkBehindIterator' that trails 'it' by ITERATOR_WINDOW_SIZE.
		static void advance_iterator_window(typename CQF<key_obj>::Iterator &it, uint64_t &step, key_obj &cqfEntry,
											typename CQF<key_obj>::Iterator &walkBehindIterator,
											const CQF<key_obj> *cqf);

		// Gathers all the distinct equivalence ID pairs and their abundance from the two DBGs into
		// the hash map 'eqClsMap'. The new IDs are set to a default value of 0.
		uint64_t gather_distinct_eq_class_pairs(ColoredDbg<qf_obj, key_obj> &dbg1, ColoredDbg<qf_obj, key_obj> &dbg2);

		// Allocates space for the buckets.
		void init_bit_vec_block_buckets(ColoredDbg<qf_obj, key_obj> &dbg1, ColoredDbg<qf_obj, key_obj> &dbg2);

		// Insert the equivalence class pair (eqCls1, eqCls2) into the hash map 'eqClsMap'.
		void add_eq_class_pair(uint64_t eqCls1, uint64_t eqCls2);

		// Concatenates the color-classes (BitVector objects) corresponding to the color-IDs 'eqID1' and 'eqID2'
		// respectively, from two CdBGs of sample sizes 'colCount1' and 'colCount2' each. The partial
		// BitVectorRRR objects containing the BitVectors for the two color-IDs are 'bv1' and 'bv2' respectively.
		// The concatenated result BitVector is stored in 'resultVec'.
		void concat(const BitVectorRRR &bv1, const uint64_t colCount1, const uint64_t eqID1,
					const BitVectorRRR &bv2, const uint64_t colCount2, const uint64_t eqID2,
					BitVector &resultVec);

		// Initializes the CQF that will contain 'finalSize' number of k-mers after mantii merge.
		void initialize_CQF(uint32_t keybits, qf_hashmode hashMode, uint32_t seed, uint64_t finalSize);

		// Concatenate required bit vectors from the DBGs into this merged DBG, to build the
		// merged color(equivalence)-class table.
		void build_eq_classes(ColoredDbg<qf_obj, key_obj> &dbg1, ColoredDbg<qf_obj, key_obj> &dbg2);

		// Merge the CQFs from the DBGs into the CQF of this DBG.
		void merge_CQFs(ColoredDbg<qf_obj, key_obj> &dbg1, ColoredDbg<qf_obj, key_obj> &dbg2);

		// Add a 'kmer' with equivalence class ID 'eqID' into the CQF of this DBG.
		void add_kmer(uint64_t kmer, uint64_t eqID, uint64_t &step,
						typename CQF<key_obj>::Iterator &walkBehindIterator);

		// Returns the sample-id mapping.
		inline std::unordered_map<uint64_t, std::string> &get_sample_id_map() { return sampleid_map; }

		// Serialize the bitvector buffer to disk, when the current class count is 'eqClsCount'.
		// TODO: Merge this overloaded method with the earlier one through default argument.
		void bv_buffer_serialize(uint64_t eqClsCount);





		// Merge approach 2 (Jamshed)

		// Name of the temporary working directory at disk; will be present temporarily inside the
		// output directory provided with the command-line arguments.
		const std::string TEMP_DIR = std::string("temp/");

		// Name of the temporary list of color-id pairs.
		const std::string EQ_ID_PAIRS_FILE = std::string("eq-id-pairs");

		// Name of the temporary file to contain count of distinct color-id pairs.
		const std::string ID_PAIR_COUNT_FILE = std::string("id-pairs-count");
		
		// Number of processor-threads to be used at the intermediate steps of 
		// unique id-pairs filtering and MPH building.
		uint threadCount = 1;

		// Maximum memory usage limit for the intermediate step of unique id-pairs filtering.
		uint maxMemoryForSort = 1;

		// Collection of distinct color-id pairs.
		std::vector<std::pair<uint64_t, uint64_t>> eqIdPair;

		// Required to hash pair objects.
		// TODO: Consult Professor on this hashing.
		// Done. Resorted to boost::hash_combine instead of plain XOR hashing.
		// For more explanation, consult https://stackoverflow.com/questions/35985960/c-why-is-boosthash-combine-the-best-way-to-combine-hash-values
		class Custom_Pair_Hasher
		{
		public:
			uint64_t operator ()(const std::pair<uint64_t, uint64_t> &key, uint64_t seed = 0) const
			{
				seed ^= std::hash<uint64_t>{}(key.first) + 0x9e3779b9 + (seed << 6) + (seed >> 2); 
				seed ^= std::hash<uint64_t>{}(key.second) + 0x9e3779b9 + (seed << 6) + (seed >> 2); 
				
				return seed;
			}
		};

		// Bloom-filter based minimal perfect hash function type.
		typedef boomphf::mphf<std::pair<uint64_t, uint64_t>, Custom_Pair_Hasher> boophf_t;

		// Minimal perfect hash function table for distinct color-id pairs.
		boophf_t *mph;

		// Redirection table for MPH(color-id-pair); required as MPH provides an arbitrary value within
		// [0, colorCount - 1] for a color-class key, whereas we need to serialize the color-ids during the
		// color-class building step, in order to reduce memory-usage by having just one output
		// bit-vector buffer for the color-class bit-vectors.
		std::vector<uint64_t> mphRedirect;

		std::vector<uint64_t> abundance;

		// Samples 'samplePairCount' number of most abundant color-id pairs from the first 'sampleKmerCount'
		// distinct k-mers of the CdBGs 'cdbg1' and 'cdbg2', into the set 'sampledPairs'.
		void sample_eq_id_pairs(ColoredDbg<qf_obj, key_obj> &cdbg1, ColoredDbg <qf_obj, key_obj> &cdbg2,
								uint64_t sampleKmerCount, uint64_t samplePairCount,
								std::unordered_set<std::pair<uint64_t, uint64_t>, Custom_Pair_Hasher> &sampledPairs);

		// Gathers all the color-id pairs for all the distinct k-mers of the CdBGs 'cdbg1' and 'cdbg2',
		// and write to a disk file; avoids repeated write of pairs present at the set 'sampledPairs'.
		uint64_t gather_eq_id_pairs(ColoredDbg<qf_obj, key_obj> &cdbg1, ColoredDbg <qf_obj, key_obj> &cdbg2,
									std::unordered_set<std::pair<uint64_t, uint64_t>, Custom_Pair_Hasher> &sampledPairs);

		// Filters all the unique color-id pairs from the pairs gathered by the method 'gather_eq_id_pairs';
		// reads and writes the id-pairs into disk.
		uint64_t gather_unique_eq_id_pairs();

		// Loads the unique id-pairs produced by the method 'gather_unique_id_pairs' into the vector 'eqIdPair'.
		void load_unique_id_pairs();

		// Sorts the vector 'eqIdPair' of unique color-id pairs based on the file-pairs they need to read from
		// in order to build the concatenated output bit-vectors.
		void block_sort();

		// Builds an MPH table 'mph' for the distinct color-id pairs present at the vector 'eqIdPair'.
		void build_mph_table();

		// Builds the output color-class bitvectors and a redirection table for the color-id pairs;
		// redirection is required as the MPH generates an arbitrary mapping of the id-pairs into the set
		// [0, colorClassCount - 1]; whereas the color-class bit-vectors are to be generated into another
		// order (the block-sorted order).
		void build_color_classes_and_redirection_table(ColoredDbg<qf_obj, key_obj> &cdbg1, ColoredDbg<qf_obj, key_obj> &cdbg2);

		// Builds the output merged CQF.
		void build_cqf(ColoredDbg<qf_obj, key_obj> &cdbg1, ColoredDbg<qf_obj, key_obj> &cdbg2);





		// Merge approach 3 (Jamshed)

		std::vector<std::vector<uint64_t>> bucketSize;

		// MPH (Minimal Perfect Hash) function tables for each of the disk-buckets.
		std::vector<std::vector<boophf_t *>> MPH;

		// The (i, j)'th entry contains the total count of color-id pairs upto bucket (i, j), exclusive,
		// in row-major order.
		std::vector<std::vector<uint64_t>> cumulativeBucketSize;

		// Initializes the disk-buckets, i.e. initializes the disk-files, MPH tables, bucket sizes,
		// cumulative size counts etc.
		inline void init_disk_buckets(ColoredDbg<qf_obj, key_obj> &cdbg1, ColoredDbg<qf_obj, key_obj> &cdbg2,
									std::vector<std::vector<std::ofstream>> &diskBucket);

		// Gathers all the color-id pairs for all the distinct k-mers of the CdBGs 'cdbg1' and 'cdbg2'
		// into disk-files (or referred to as buckets hereafter), where an id-pair goes to bucket_(i, j)
		// if it reads color-class bitvectors from bitvector_file_(i-1) of cdbg1 and bitvector_file_(j-1)
		// of cdbg2; with avoiding repeated writes of pairs from the set 'sampledPairs', sampled on abundance.
		// A bucket of the form (0, X) with X > 0 implies that, the id-pairs present at this bucket are only
		// of k-mers that are absent at cdbg1, present at cdbg2, and read from the bitvector_file_(X - 1)
		// of cdbg2. Buckets of the form (X, 0) imply vice versa.
		// Returns the number of distinct k-mers present at the CdBGs cdg1 and cdbg2.
		uint64_t fill_disk_buckets(ColoredDbg<qf_obj, key_obj> &cdbg1, ColoredDbg<qf_obj, key_obj> &cdbg2,
								std::unordered_set<std::pair<uint64_t, uint64_t>, Custom_Pair_Hasher> &sampledPairs);

		// Adds the color-class ID pair (colorID1, colorID2) to the appropriate disk bucket;
		// i.e. write the pair into the file diskBucket[i][j] iff colorID1 has its color-class (bitvector)
		// at bitvector_file_i and colorID2 has its color-class at bitvector_file_j.
		inline void add_color_id_pair(const uint64_t colorID1, const uint64_t colorID2,
										std::vector<std::vector<std::ofstream>> &diskBucket);

		// Filters all the disk-buckets filled up by the method 'fill_disk_buckets', to contain only unique
		// color-id pairs. Returns the count of unique color-id pairs.
		uint64_t filter_disk_buckets(ColoredDbg<qf_obj, key_obj> &cdbg1, ColoredDbg<qf_obj, key_obj> &cdbg2);

		// Builds an MPH (Minimal Perfect Hash) table for each disk-bucket.
		void build_MPH_tables(ColoredDbg<qf_obj, key_obj> &cdbg1, ColoredDbg<qf_obj, key_obj> &cdbg2);

		// Builds the output color-class bitvectors for the color-id pairs.
		void build_color_class_table(ColoredDbg<qf_obj, key_obj> &cdbg1, ColoredDbg<qf_obj, key_obj> &cdbg2);

		// Given a color-id pair 'idPair', returns the newly assigned color-id of this pair,
		// at the merged CdBG.
		inline uint64_t get_color_id(const std::pair<uint64_t, uint64_t> &idPair);

		// Builds the output merged CQF.
		void build_CQF(ColoredDbg<qf_obj, key_obj> &cdbg1, ColoredDbg<qf_obj, key_obj> &cdbg2);
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
      walk_behind_iterator = dbg.begin(true);
    } else if (counter > 4096) {
      ++walk_behind_iterator;
    }
    
		// Progress tracker
		static uint64_t last_size = 0;
		if (dbg.dist_elts() % PROGRESS_STEP == 0 &&
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
ColoredDbg<qf_obj, key_obj>::
	ColoredDbg(ColoredDbg<qf_obj, key_obj> &cdbg1, ColoredDbg<qf_obj, key_obj> &cdbg2,
				std::string &prefix, int flag):
	bv_buffer(), prefix(prefix), num_samples(cdbg1.get_num_samples() + cdbg2.get_num_samples()),
	num_serializations(0), threadCount(1), dbg_alloc_flag(flag),
	start_time_(std::time(nullptr))
	{
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

		while(sampleList >> sampleID >> sampleName)
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

			// printf("\nLoaded color-class file %s\n", file.c_str());
		}
		

		eqClsFiles.reserve(eqclassFiles.size());
		for(auto idFilePair : sortedFiles)
			eqClsFiles.push_back(idFilePair.second);

		num_serializations = eqClsFiles.size();

		// puts("MSG: Done loading colored dBG.\n");
}





template <typename qf_obj, typename key_obj>
inline void ColoredDbg<qf_obj, key_obj> ::
	init_bit_vec_block_buckets(ColoredDbg<qf_obj, key_obj> &dbg1, ColoredDbg<qf_obj, key_obj> &dbg2)
{
	const uint64_t fileCount1 = dbg1.get_eq_class_file_count(), fileCount2 = dbg2.get_eq_class_file_count();
	
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

		uint64_t file1 = (eqCls1 ? (eqCls1 - 1) / mantis::NUM_BV_BUFFER : 0),
					file2 = (eqCls2 ? (eqCls2 - 1) / mantis::NUM_BV_BUFFER : 0);

		bucket[file1][file2].push_back(eqPair);
	}
	else
		it -> second.second++;
}





template <typename qf_obj, typename key_obj>
inline void ColoredDbg<qf_obj, key_obj> ::
	advance_iterator_window(typename CQF<key_obj>::Iterator &it, uint64_t &step, key_obj &cqfEntry,
							typename CQF<key_obj>::Iterator &walkBehindIterator, const CQF<key_obj> *cqf)
{
	++it, ++step;
	if(!it.done())
		cqfEntry = it.get_cur_hash();

	if(step == ITERATOR_WINDOW_SIZE)
		walkBehindIterator = cqf -> begin(true);
	else if(step > ITERATOR_WINDOW_SIZE)
		++walkBehindIterator;
}





template <typename qf_obj, typename key_obj>
uint64_t ColoredDbg<qf_obj, key_obj> ::
	gather_distinct_eq_class_pairs(ColoredDbg<qf_obj, key_obj> &dbg1, ColoredDbg<qf_obj, key_obj> &dbg2)
{
	auto t_start = time(nullptr);


	console -> info("At distinct equivalence-class pairs gathering phase. Time = {}", time(nullptr) - start_time_);

	const CQF<key_obj> *cqf1 = dbg1.get_cqf(), *cqf2 = dbg2.get_cqf();
	typename CQF<key_obj>::Iterator it1 = cqf1 -> begin(), it2 = cqf2 -> begin(),
									walkBehindIterator1, walkBehindIterator2;
	uint64_t cqfPosition1 = 0, cqfPosition2 = 0;

	uint64_t kmerCount = 0, equalKmerCount = 0;
	uint64_t kmer1, kmer2, eqClass1, eqClass2;
	key_obj cqfEntry1, cqfEntry2;

	
	init_bit_vec_block_buckets(dbg1, dbg2);


	if(it1.done() || it2.done())
		console -> error("One or more CQF iterator(s) already at end position before starting walk.");


	if(!it1.done())
		cqfEntry1 = it1.get_cur_hash(), cqfPosition1++;
	
	if(!it2.done())
		cqfEntry2 = it2.get_cur_hash(), cqfPosition2++;


	while(!it1.done() || !it2.done())
	{
		if(!it1.done())
			kmer1 = cqfEntry1.key;
		
		if(!it2.done())
			kmer2 = cqfEntry2.key;


		// eqClassX = 0 implies absence in CdBG X.
		if(it1.done())
		{
			eqClass1 = 0, eqClass2 = cqfEntry2.count;
			advance_iterator_window(it2, cqfPosition2, cqfEntry2, walkBehindIterator2, cqf2);// ++it2;
		}
		else if(it2.done())
		{
			eqClass1 = cqfEntry1.count, eqClass2 = 0;
			advance_iterator_window(it1, cqfPosition1, cqfEntry1, walkBehindIterator1, cqf1);// ++it1;
		}
		else if(kmer1 < kmer2)
		{
			eqClass1 = cqfEntry1.count, eqClass2 = 0;
			advance_iterator_window(it1, cqfPosition1, cqfEntry1, walkBehindIterator1, cqf1);// ++it1;
		}
		else if(kmer2 < kmer1)
		{
			eqClass1 = 0, eqClass2 = cqfEntry2.count;
			advance_iterator_window(it2, cqfPosition2, cqfEntry2, walkBehindIterator2, cqf2);// ++it2;
		}
		else
		{
			eqClass1 = cqfEntry1.count, eqClass2 = cqfEntry2.count;
			advance_iterator_window(it1, cqfPosition1, cqfEntry1, walkBehindIterator1, cqf1),// ++it1;
			advance_iterator_window(it2, cqfPosition2, cqfEntry2, walkBehindIterator2, cqf2),// ++it2;
			equalKmerCount++;
		}


		add_eq_class_pair(eqClass1, eqClass2);
		kmerCount++;


		if(kmerCount % PROGRESS_STEP == 0)
			console -> info("Observed count of -- distinct k-mers: {}M, shared k-mers: {}M, color-class: {}. Time {}",
							kmerCount * 10 / PROGRESS_STEP, equalKmerCount * 10.0 / PROGRESS_STEP,
							eqClsMap.size(), time(nullptr) - start_time_);
	}


	console -> info("Distinct kmers found {}, shared kmers found {}, color-class count {}. Time = {}",
					kmerCount, equalKmerCount, eqClsMap.size(), time(nullptr) - start_time_);

	
	auto t_end = time(nullptr);
	console -> info("Phase 1 took time {} seconds.", t_end - t_start);

	return kmerCount;
}





template <typename qf_obj, typename key_obj>
void ColoredDbg<qf_obj, key_obj> ::
	initialize_CQF(uint32_t keybits, qf_hashmode hashMode, uint32_t seed, uint64_t finalSize)
{
	// Get floor(log2(finalSize))
	uint32_t qbits;
	for(qbits = 0; (finalSize >> qbits) != (uint64_t)1; qbits++);

	// Get ceil(log2(finalSize))
	if(finalSize & (finalSize - 1))	// if finalSize is not a power of 2
		qbits++;

	
	qbits += 3;	// to avoid the initial rapid resizes at minuscule load factors

	
	if(dbg_alloc_flag == MANTIS_DBG_IN_MEMORY)
	{
		CQF<key_obj> cqf(qbits, keybits, hashMode, seed);
		dbg = cqf;
	}
	else if(dbg_alloc_flag == MANTIS_DBG_ON_DISK)
	{
		CQF<key_obj> cqf(qbits, keybits, hashMode, seed, prefix + mantis::CQF_FILE);
		dbg = cqf;
	}
	else
	{
		ERROR("Wrong Mantis alloc mode.");
		exit(EXIT_FAILURE);
	}

	dbg.set_auto_resize();
}





template <typename qf_obj, typename key_obj>
void ColoredDbg<qf_obj, key_obj> ::
	concat(const BitVectorRRR &bv1, const uint64_t colCount1, const uint64_t eqID1,
			const BitVectorRRR &bv2, const uint64_t colCount2, const uint64_t eqID2,
			BitVector &resultVec)
{
	const uint64_t wordLen = 64;


	if(eqID1)	// Eq ID = 0 implies an absent equivalence-id.
	{
		// Optimized bitvector read and write

		uint64_t offset = ((eqID1 - 1) % mantis::NUM_BV_BUFFER) * colCount1;

		for(uint32_t blockStart = 0; blockStart < (colCount1 / wordLen) * wordLen; blockStart += wordLen)
				resultVec.set_int(blockStart, bv1.get_int(offset + blockStart, wordLen), wordLen);

		if(colCount1 % wordLen)
			resultVec.set_int((colCount1 / wordLen) * wordLen,
								bv1.get_int(offset + ((colCount1 / wordLen) * wordLen), colCount1 % wordLen),
								colCount1 % wordLen);



		// Slower version (should be, uses [] operator)

		// bool isEmpty = true; // for debugging purpose(s)

		// uint64_t offset = ((eqID1 - 1) % mantis::NUM_BV_BUFFER) * colCount1;

		// for(uint32_t wordCount = 0; wordCount <= colCount1 / wordLen; ++wordCount)
		// {
		// 	uint64_t readLen = std::min(wordLen, colCount1 - wordCount * wordLen);
		// 	uint64_t word = bv1.get_int(offset, readLen);

		// 	// printf("eqId1 %d, read length = %d, word %d\n", (int)eqID1, (int)readLen, (int)word);

		// 	// Optimize here; preferrably eliminate the loop with one statement (some sort of set_int() ?).
		// 	for(uint32_t bitIdx = 0, sampleID = wordCount * wordLen; bitIdx < readLen; ++bitIdx, ++sampleID)
		// 		if((word >> bitIdx) & 0x01)
		// 			resultVec[sampleID] = 1, isEmpty = false;			

		// 	offset += readLen;
		// }
	}


	if(eqID2)	// Eq ID = 0 implies an absent equivalence id.
	{
		// Optimized bitvector read and write

		uint64_t offset = ((eqID2 - 1) % mantis::NUM_BV_BUFFER) * colCount2;

		for(uint32_t blockStart = 0; blockStart < (colCount2 / wordLen) * wordLen; blockStart += wordLen)
				resultVec.set_int(colCount1 + blockStart, bv2.get_int(offset + blockStart, wordLen), wordLen);

		if(colCount2 % wordLen)
			resultVec.set_int(colCount1 + (colCount2 / wordLen) * wordLen,
								bv2.get_int(offset + ((colCount2 / wordLen) * wordLen), colCount2 % wordLen),
								colCount2 % wordLen);


		// Slower version (should be, uses [] operator)

		// u_int64_t offset2 = ((eqID2 - 1) % mantis::NUM_BV_BUFFER) * colCount2;

		// for(uint32_t wordCount = 0; wordCount <= colCount2 / wordLen; ++wordCount)
		// {
		// 	uint64_t readLen = std::min(wordLen, colCount2 - wordCount * wordLen);
		// 	uint64_t word = bv2.get_int(offset2, readLen);

		// 	// printf("eqId2 %d, read length = %d, word %d\n", (int)eqID2, (int)readLen, (int)word);

		// 	// Optimize here; preferrably eliminate the loop with one statement (some sort of set_int() ?).
		// 	for(uint32_t bitIdx = 0, sampleID = wordCount * wordLen; bitIdx < readLen; ++bitIdx, ++sampleID)
		// 		if((word >> bitIdx) & 0x01)
		// 			resultVec[colCount1 + sampleID] = 1;

		// 	offset2 += readLen;
		// }
	}
}





template <class qf_obj, class key_obj>
void ColoredDbg<qf_obj, key_obj>::bv_buffer_serialize(uint64_t eqClsCount)
{
	// BitVector bv_temp(bv_buffer);
	if (eqClsCount % mantis::NUM_BV_BUFFER > 0)
		bv_buffer.resize((eqClsCount % mantis::NUM_BV_BUFFER) * num_samples);

	
	BitVectorRRR final_com_bv(bv_buffer);
	std::string bv_file(prefix + std::to_string(num_serializations) + "_" + mantis::EQCLASS_FILE);
	

	sdsl::store_to_file(final_com_bv, bv_file);
	// bv_buffer = BitVector(bv_buffer.bit_size());	// why is it required?

	num_serializations++;
}





template <typename qf_obj, typename key_obj>
void ColoredDbg<qf_obj, key_obj> ::
	build_eq_classes(ColoredDbg<qf_obj, key_obj> &dbg1, ColoredDbg<qf_obj, key_obj> &dbg2)
{
	auto t_start = time(nullptr);


	console -> info("At color-class building (bitvectors concatenation) phase. Time = {}\n",
					time(nullptr) - start_time_);

	const uint64_t fileCount1 = dbg1.get_eq_class_file_count(),
					fileCount2 = dbg2.get_eq_class_file_count();
	uint64_t serialID = 0;
	


	for(uint64_t i = 0; i < fileCount1; ++i)
	{
		console -> info("Mantis 1: loading one bitvectorRRR from file {}.", dbg1.get_eq_class_files()[i]);

		// Required: data_1 = read_i'th Bit Vector block for dbg1
		BitVectorRRR bitVec1;
		sdsl::load_from_file(bitVec1, dbg1.get_eq_class_files()[i]);

		// printf("\nLoading bitvectorRRR from file %s.\n", dbg1.get_eq_class_files()[i].c_str());

		for(uint64_t j = 0; j < fileCount2; ++j)
		{
			console -> info("At bucket ({}, {}). Size = {}", i, j, bucket[i][j].size());

			console -> info("Mantis 2: loading one bitvectorRRR from file {}.", dbg2.get_eq_class_files()[j]);
			

			// Required: data_2 = read_j'th Bit Vector block for dbg2
			BitVectorRRR bitVec2;
			sdsl::load_from_file(bitVec2, dbg2.get_eq_class_files()[j]);

			// printf("\nLoading bitvectorRRR from file %s\n", dbg2.get_eq_class_files()[j].c_str());

			for(auto eqClsPair : bucket[i][j])
			{
				uint64_t eq1 = eqClsPair.first, eq2 = eqClsPair.second;
				eqClsMap[eqClsPair].first = ++serialID; // serialID++ doesn't work. Introduces bug in CQF;
														// Probable cause: CQF doesn't support 0-count ?

				// Required: bit_vector = concat(data_1.query(eq1), data_2.query(eq2))
				// Required: dbg.bitVectorBuffer[serialID] = bit_vector
				BitVector mergedEqCls(num_samples);
				concat(bitVec1, dbg1.get_num_samples(), eq1, bitVec2, dbg2.get_num_samples(), eq2, mergedEqCls);

				
				// for debugging purpose(s)
				// printf("bitvector for eq id %d is = ", (int)serialID);
				// bool isEmpty = true;
				// for(int i = 0; i < num_samples; ++i)
				// 	if(mergedEqCls[i] == 1)
				// 		putchar('1'), isEmpty = false;
				// 	else
				// 		putchar('0');
				// putchar('\n');

				// puts("Going to add a bitvector\n");
				add_bitvector(mergedEqCls, serialID - 1);
				// puts("Bitvector added\n");

				// if(isEmpty)
				// 	puts("\nEmpty bitvector constructed\n");
				// else
				// 	puts("Nonempty color class.");

				
				// Serialization and disk-write if required
				if(serialID % mantis::NUM_BV_BUFFER == 0)
				{
					console -> info("Serializing bitvector buffer with {} color-classes.", serialID);
					bv_buffer_serialize(serialID);
				}
			}


			// Clear memory for this bucket.
			bucket[i][j].clear();
			bucket[i][j].shrink_to_fit();
		}
	}


	// Serialize the bv buffer last time if needed
	if (serialID % mantis::NUM_BV_BUFFER > 0)
	{
		console -> info("Serializing bitvector buffer with {} color-classes.", serialID);
		bv_buffer_serialize(serialID);
	}

	
	auto t_end = time(nullptr);
	console -> info("Phase 2 took time {} seconds.", t_end - t_start);
}





template <typename qf_obj, typename key_obj>
inline void ColoredDbg<qf_obj, key_obj> ::
	add_kmer(uint64_t kmer, uint64_t eqID, uint64_t &step, typename CQF<key_obj>::Iterator &walkBehindIterator)
{
	if(dbg.insert(KeyObject(kmer, 0, eqID), QF_NO_LOCK | QF_KEY_IS_HASH) == QF_NO_SPACE)
	{
		// This means that auto_resize failed.
		console -> error("The CQF is full and auto resize failed. Please rerun build with a bigger size.");
		exit(1);
	}

	
	step++;
	if(step == ITERATOR_WINDOW_SIZE)
		walkBehindIterator = dbg.begin(true);
	else if(step > ITERATOR_WINDOW_SIZE)
		++walkBehindIterator;
}





template <typename qf_obj, typename key_obj>
void ColoredDbg<qf_obj, key_obj> ::
	merge_CQFs(ColoredDbg<qf_obj, key_obj> &dbg1, ColoredDbg<qf_obj, key_obj> &dbg2)
{
	auto t_start = time(nullptr);


	console -> info("At CQFs merging phase. Time = {}\n", time(nullptr) - start_time_);

	const CQF<key_obj> *cqf1 = dbg1.get_cqf(), *cqf2 = dbg2.get_cqf();
	typename CQF<key_obj>::Iterator it1 = cqf1 -> begin(), it2 = cqf2 -> begin(),
									walkBehindIterator1, walkBehindIterator2, walkBehindIteratorOut;
	uint64_t cqfPosition1 = 0, cqfPosition2 = 0, cqfOutPosition = 0;

	uint64_t kmerCount = 0, equalKmerCount = 0;
	uint64_t kmer1, kmer2, kmer, eqClass1, eqClass2;
	key_obj cqfEntry1, cqfEntry2;

	
	if(it1.done() || it2.done())
		console -> error("One or more CQF iterator(s) already at end position before starting walk.");

	
	if(!it1.done())
		cqfEntry1 = it1.get_cur_hash(), cqfPosition1++;
	
	if(!it2.done())
		cqfEntry2 = it2.get_cur_hash(), cqfPosition2++;


	while(!it1.done() || !it2.done())
	{
		if(!it1.done())
			kmer1 = cqfEntry1.key;
		
		if(!it2.done())
			kmer2 = cqfEntry2.key;


		if(it1.done())
		{
			kmer = kmer2;

			eqClass1 = 0, eqClass2 = cqfEntry2.count;
			advance_iterator_window(it2, cqfPosition2, cqfEntry2, walkBehindIterator2, cqf2);// ++it2;
		}
		else if(it2.done())
		{
			kmer = kmer1;
			
			eqClass1 = cqfEntry1.count, eqClass2 = 0;
			advance_iterator_window(it1, cqfPosition1, cqfEntry1, walkBehindIterator1, cqf1);// ++it1;
		}
		else if(kmer1 < kmer2)
		{
			kmer = kmer1;
			
			eqClass1 = cqfEntry1.count, eqClass2 = 0;
			advance_iterator_window(it1, cqfPosition1, cqfEntry1, walkBehindIterator1, cqf1);// ++it1;
		}
		else if(kmer2 < kmer1)
		{
			kmer = kmer2;

			eqClass1 = 0, eqClass2 = cqfEntry2.count;
			advance_iterator_window(it2, cqfPosition2, cqfEntry2, walkBehindIterator2, cqf2);// ++it2;
		}
		else
		{
			kmer = kmer1;

			eqClass1 = cqfEntry1.count, eqClass2 = cqfEntry2.count;
			advance_iterator_window(it1, cqfPosition1, cqfEntry1, walkBehindIterator1, cqf1),// ++it1;
			advance_iterator_window(it2, cqfPosition2, cqfEntry2, walkBehindIterator2, cqf2);// ++it2;

			equalKmerCount++; // for debugging purpose(s)
		}


		add_kmer(kmer, eqClsMap[std :: make_pair(eqClass1, eqClass2)].first,
					cqfOutPosition, walkBehindIteratorOut);
		kmerCount++;

		if(kmerCount % PROGRESS_STEP == 0)
			console -> info("Kmers merged: {}M, time: {}", kmerCount  * 10 / PROGRESS_STEP,
							time(nullptr) - start_time_);
	}


	console -> info("Total kmers merged: {}. Time: {}", kmerCount, time(nullptr) - start_time_);


	auto t_end = time(nullptr);
	console -> info("Phase 3 took time {} seconds.", t_end - t_start);
}





template <typename qf_obj, typename key_obj>
void ColoredDbg<qf_obj, key_obj> :: 
	construct(ColoredDbg<qf_obj, key_obj> &dbg1, ColoredDbg<qf_obj, key_obj> &dbg2)
{
	console -> info ("Merge starting. Time = {}\n", time(nullptr) - start_time_);
	

	uint64_t kmerCount = gather_distinct_eq_class_pairs(dbg1, dbg2);

	build_eq_classes(dbg1, dbg2);

	initialize_CQF(dbg1.get_cqf() -> keybits(), dbg1.get_cqf() -> hash_mode(), dbg1.get_cqf() -> seed(),
					kmerCount);
	merge_CQFs(dbg1, dbg2);
	

	console -> info ("Merge ending. Time = {}\n", time(nullptr) - start_time_);
}





template <typename qf_obj, typename key_obj>
void ColoredDbg<qf_obj, key_obj>::
	serialize(ColoredDbg<qf_obj, key_obj> &dbg1, ColoredDbg<qf_obj, key_obj> &dbg2)
{
	// Serialize the CQF
	if (dbg_alloc_flag == MANTIS_DBG_IN_MEMORY)
		dbg.serialize(prefix + mantis::CQF_FILE);
	else
		dbg.close();


	// Serialize the bv buffer last time if needed.
	// Done at the color(equivalence)-class building phase.


	// TODO: Move the list concatenation phase to constructor, with introduction of a new field;
	// just keep the disk-write part here.

	// Serialize the sample id map
	std::ofstream outputFile(prefix + mantis::SAMPLEID_FILE);

	for(auto sampleID: dbg1.get_sample_id_map())
		outputFile << sampleID.first << " " << sampleID.second << "\n";

	for(auto sampleID: dbg2.get_sample_id_map())
		outputFile << dbg1.get_num_samples() + sampleID.first << " " << sampleID.second << "\n";

	outputFile.close();


	// Dump the abundance distribution of the equivalence / color classes.
	if (flush_eqclass_dis)
	{
		const char OUTPUT_ABUNDANCE_DIST_FILE[] = "eqclass_dist.lst";
		std::ofstream outputFile(prefix + OUTPUT_ABUNDANCE_DIST_FILE);

		for (auto idFreq : eqClsMap)
			outputFile << idFreq.second.first << " " << idFreq.second.second << "\n";
		
		outputFile.close();
	}
}





template <typename qf_obj, typename key_obj>
void ColoredDbg<qf_obj, key_obj>::
	sample_eq_id_pairs(ColoredDbg<qf_obj, key_obj> &cdbg1, ColoredDbg <qf_obj, key_obj> &cdbg2,
	uint64_t sampleKmerCount, uint64_t samplePairCount,
	std::unordered_set<std::pair<uint64_t, uint64_t>, Custom_Pair_Hasher> &sampledPairs)
{
	auto t_start = time(nullptr);

	console -> info("Sampling {} most-abundant eq-id pairs from first {} kmers. Time-stamp = {}",
					samplePairCount, sampleKmerCount, time(nullptr) - start_time_);

	
	const CQF<key_obj> *cqf1 = cdbg1.get_cqf(), *cqf2 = cdbg2.get_cqf();
	typename CQF<key_obj>::Iterator it1 = cqf1 -> begin(), it2 = cqf2 -> begin(),
									walkBehindIterator1, walkBehindIterator2;
	uint64_t cqfPosition1 = 0, cqfPosition2 = 0;

	uint64_t kmerCount = 0;
	uint64_t kmer1, kmer2, eqClass1, eqClass2;
	key_obj cqfEntry1, cqfEntry2;

	std::unordered_map<std::pair<uint64_t, uint64_t>, uint64_t, Custom_Pair_Hasher> pairCount;


	if(it1.done() || it2.done())
		console -> error("One or more CQF iterator(s) already at end position before starting walk.");


	if(!it1.done())
		cqfEntry1 = it1.get_cur_hash(), cqfPosition1++;
	
	if(!it2.done())
		cqfEntry2 = it2.get_cur_hash(), cqfPosition2++;


	while(kmerCount < sampleKmerCount && (!it1.done() || !it2.done()))
	{
		if(!it1.done())
			kmer1 = cqfEntry1.key;
		
		if(!it2.done())
			kmer2 = cqfEntry2.key;


		// eqClassX = 0 implies absence in CdBG X.
		if(it1.done())
		{
			eqClass1 = 0, eqClass2 = cqfEntry2.count;
			advance_iterator_window(it2, cqfPosition2, cqfEntry2, walkBehindIterator2, cqf2);// ++it2;
		}
		else if(it2.done())
		{
			eqClass1 = cqfEntry1.count, eqClass2 = 0;
			advance_iterator_window(it1, cqfPosition1, cqfEntry1, walkBehindIterator1, cqf1);// ++it1;
		}
		else if(kmer1 < kmer2)
		{
			eqClass1 = cqfEntry1.count, eqClass2 = 0;
			advance_iterator_window(it1, cqfPosition1, cqfEntry1, walkBehindIterator1, cqf1);// ++it1;
		}
		else if(kmer2 < kmer1)
		{
			eqClass1 = 0, eqClass2 = cqfEntry2.count;
			advance_iterator_window(it2, cqfPosition2, cqfEntry2, walkBehindIterator2, cqf2);// ++it2;
		}
		else
		{
			eqClass1 = cqfEntry1.count, eqClass2 = cqfEntry2.count;
			advance_iterator_window(it1, cqfPosition1, cqfEntry1, walkBehindIterator1, cqf1),// ++it1;
			advance_iterator_window(it2, cqfPosition2, cqfEntry2, walkBehindIterator2, cqf2);// ++it2;
		}


		pairCount[std::make_pair(eqClass1, eqClass2)]++;
		kmerCount++;


		if(kmerCount % PROGRESS_STEP == 0)
			console -> info("Sampled {}M k-mers, color-classes found: {}. Time-stamp = {}",
							kmerCount * 10 / PROGRESS_STEP, pairCount.size(),
							time(nullptr) - start_time_);
	}

	
	console -> info("Sampled {} k-mers, color-classes found: {}. Time-stamp = {}",
					kmerCount, pairCount.size(), time(nullptr) - start_time_);



	typedef std::pair<uint64_t, std::pair<uint64_t, uint64_t>> CountAndIdPair;
	std::priority_queue<CountAndIdPair, std::vector<CountAndIdPair>, std::greater<CountAndIdPair>> minPQ;

	for(auto p = pairCount.begin(); p != pairCount.end(); ++p)
		if(minPQ.size() < samplePairCount)
			minPQ.push(std::make_pair(p -> second, p -> first));
		else if(minPQ.top().first < p -> second)
		{
			minPQ.pop();
			minPQ.push(std::make_pair(p -> second, p -> first));
		}


	
	while(!minPQ.empty())
	{
		sampledPairs.insert(minPQ.top().second);
		minPQ.pop();
	}


	console -> info("Sampled {} eq-id pairs. Time-stamp = {}", sampledPairs.size(), time(nullptr) - start_time_);


	auto t_end = time(nullptr);
	console -> info("Sampling abundant equivalence-id pairs took time {} seconds.", t_end - t_start);
}






template <typename qf_obj, typename key_obj>
uint64_t ColoredDbg<qf_obj, key_obj>::
	gather_eq_id_pairs(ColoredDbg<qf_obj, key_obj> &cdbg1, ColoredDbg <qf_obj, key_obj> &cdbg2,
						std::unordered_set<std::pair<uint64_t, uint64_t>, Custom_Pair_Hasher> &sampledPairs)
{
	auto t_start = time(nullptr);


	console -> info("Writing all the eq-id pairs to file {}. Time-stamp = {}", TEMP_DIR + EQ_ID_PAIRS_FILE,
					time(nullptr) - start_time_);


	
	// TODO: Add faster file-write mechanism.
	std::string fileName = prefix + TEMP_DIR + EQ_ID_PAIRS_FILE;
	std::ofstream output(fileName);


	uint64_t writtenPairsCount;

	for(auto p = sampledPairs.begin(); p != sampledPairs.end(); ++p)
		output << p -> first << " " << p -> second << "\n";

	writtenPairsCount = sampledPairs.size();

	console -> info("Wrote {} sampled pairs. Time-stamp = {}", writtenPairsCount, time(nullptr) - start_time_);


	console -> info("Now iterating over the CQFs for the rest of the id pairs.");

	const CQF<key_obj> *cqf1 = cdbg1.get_cqf(), *cqf2 = cdbg2.get_cqf();
	typename CQF<key_obj>::Iterator it1 = cqf1 -> begin(), it2 = cqf2 -> begin(),
									walkBehindIterator1, walkBehindIterator2;
	uint64_t cqfPosition1 = 0, cqfPosition2 = 0;

	uint64_t kmerCount = 0;
	uint64_t kmer1, kmer2, eqClass1, eqClass2;
	key_obj cqfEntry1, cqfEntry2;


	if(it1.done() || it2.done())
		console -> error("One or more CQF iterator(s) already at end position before starting walk.");


	if(!it1.done())
		cqfEntry1 = it1.get_cur_hash(), cqfPosition1++;
	
	if(!it2.done())
		cqfEntry2 = it2.get_cur_hash(), cqfPosition2++;


	while(!it1.done() || !it2.done())
	{
		if(!it1.done())
			kmer1 = cqfEntry1.key;
		
		if(!it2.done())
			kmer2 = cqfEntry2.key;


		// eqClassX = 0 implies absence in CdBG X.
		if(it1.done())
		{
			eqClass1 = 0, eqClass2 = cqfEntry2.count;
			advance_iterator_window(it2, cqfPosition2, cqfEntry2, walkBehindIterator2, cqf2);// ++it2;
		}
		else if(it2.done())
		{
			eqClass1 = cqfEntry1.count, eqClass2 = 0;
			advance_iterator_window(it1, cqfPosition1, cqfEntry1, walkBehindIterator1, cqf1);// ++it1;
		}
		else if(kmer1 < kmer2)
		{
			eqClass1 = cqfEntry1.count, eqClass2 = 0;
			advance_iterator_window(it1, cqfPosition1, cqfEntry1, walkBehindIterator1, cqf1);// ++it1;
		}
		else if(kmer2 < kmer1)
		{
			eqClass1 = 0, eqClass2 = cqfEntry2.count;
			advance_iterator_window(it2, cqfPosition2, cqfEntry2, walkBehindIterator2, cqf2);// ++it2;
		}
		else
		{
			eqClass1 = cqfEntry1.count, eqClass2 = cqfEntry2.count;
			advance_iterator_window(it1, cqfPosition1, cqfEntry1, walkBehindIterator1, cqf1),// ++it1;
			advance_iterator_window(it2, cqfPosition2, cqfEntry2, walkBehindIterator2, cqf2);// ++it2;
		}


		kmerCount++;

		if(sampledPairs.find(std::make_pair(eqClass1, eqClass2)) == sampledPairs.end())
		{
			output << eqClass1 << " " << eqClass2 << "\n";
			writtenPairsCount++;
		}


		if(kmerCount % PROGRESS_STEP == 0)
			console -> info("Observed count of distinct k-mers: {}M, written id-pairs: {}. Time-stamp = {}",
							kmerCount * 10 / PROGRESS_STEP, writtenPairsCount, time(nullptr) - start_time_);
	}


	console -> info("Distinct kmers found {}, id pairs written {}. Time-stamp = {}",
					kmerCount, writtenPairsCount, time(nullptr) - start_time_);

	output.flush();
	output.close();

	
	auto t_end = time(nullptr);
	console -> info("Gathering equivalence-id pairs took time {} seconds.", t_end - t_start);


	return kmerCount;
}





template <typename qf_obj, typename key_obj>
uint64_t ColoredDbg<qf_obj, key_obj>:: gather_unique_eq_id_pairs()
{
	if(!system(NULL))
	{
		console -> error("Command processor does not exist.");
		exit(1);
	}


	auto t_start = time(nullptr);

	console -> info("Filtering out the unique eq-id pairs from file {} with {} threads. Time-stamp = {}",
					TEMP_DIR + EQ_ID_PAIRS_FILE, threadCount, time(nullptr) - start_time_);


	std::string pairsFile = prefix + TEMP_DIR + EQ_ID_PAIRS_FILE;
	std::string sysCommand = "sort -u";

	// TODO: Consult professor on parallelization and memory-usage details.
	sysCommand += " --parallel=" + std::to_string(threadCount);
	sysCommand += " -S " + std::to_string(maxMemoryForSort) + "G";
	sysCommand += " -o " + pairsFile + " " + pairsFile;

	console -> info("System command used:\n{}", sysCommand);

	system(sysCommand.c_str());

	console -> info("Filtered all the unique eq-id pairs into file {}. Time-stamp = {}",
					pairsFile, time(nullptr) - start_time_);



	std::string countFile = prefix + TEMP_DIR + ID_PAIR_COUNT_FILE;
	sysCommand = "wc -l " + pairsFile + " | egrep -o \"[0-9]+ \" > " + countFile;

	system(sysCommand.c_str());

	std::ifstream input(countFile);
	uint64_t colorClassCount;

	input >> colorClassCount;
	input.close();

	console -> info("Count of unique pairs = {}. Time-stamp = {}", colorClassCount, time(nullptr) - start_time_);

	
	auto t_end = time(nullptr);
	console -> info("Filtering the unique id pairs took time {} seconds.", t_end - t_start);

	return colorClassCount;
}





template <typename qf_obj, typename key_obj>
inline void ColoredDbg<qf_obj, key_obj>:: load_unique_id_pairs()
{
	auto t_start = time(nullptr);

	console -> info("Loading the unique id pairs into memory.");


	std::string pairsFile = prefix + TEMP_DIR + EQ_ID_PAIRS_FILE;
	std::ifstream input(pairsFile);
	uint64_t id1, id2;;

	while(input >> id1 >> id2)
		eqIdPair.emplace_back(id1, id2);

	input.close();

	console -> info("Read {} id pairs into memory. Time-stamp = {}", eqIdPair.size(),
					time(nullptr) - start_time_);


	auto t_end = time(nullptr);
	console -> info("Loading the id pairs into memory took time {} seconds.", t_end - t_start);
}





template <typename qf_obj, typename key_obj>
inline void ColoredDbg<qf_obj, key_obj>:: block_sort()
{
	auto t_start = time(nullptr);


	console -> info("Sorting {} unique eq-id pairs based on their corresponding color-class blocks. Time-stamp = {}",
					eqIdPair.size(), time(nullptr) - start_time_);

	sort(eqIdPair.begin(), eqIdPair.end(),
		[](const std::pair<uint64_t, uint64_t> &lhs, const std::pair<uint64_t, uint64_t> &rhs)
		{
			uint64_t bucket1 = (lhs.first ? (lhs.first - 1) / mantis::NUM_BV_BUFFER : 0),
						bucket2 = (rhs.first ? (rhs.first - 1) / mantis::NUM_BV_BUFFER : 0);

			if(bucket1 != bucket2)
				return bucket1 < bucket2;

			bucket1 = (lhs.second ? (lhs.second - 1) / mantis::NUM_BV_BUFFER : 0),
			bucket2 = (rhs.second ? (rhs.second - 1) / mantis::NUM_BV_BUFFER : 0);

			return bucket1 < bucket2;
		});

	console -> info("Sorting complete. Time-stamp = {}", time(nullptr) - start_time_);


	auto t_end = time(nullptr);
	console -> info("Block-sorting the unique equivalence-id pairs took time {} seconds.", t_end - t_start);
}





template <typename qf_obj, typename key_obj>
void ColoredDbg<qf_obj, key_obj>:: build_mph_table()
{
	auto t_start = time(nullptr);

	console -> info("Building the minimal perfect hashing table. Time-stamp = {}", time(nullptr) - start_time_);
	

	uint64_t nelem = eqIdPair.size();
	uint nthreads = threadCount;
	
	mph = NULL;
	console -> info("Constructing a BooPHF with {} elements using {} threads.", nelem, nthreads);
	
	
	// lowest bit/elem is achieved with gamma=1, higher values lead to larger mphf but faster construction/query
	double gammaFactor = 2.0;	// gamma = 2 is a good tradeoff (leads to approx 3.7 bits/key )

	// build the mphf
	mph = new boomphf::mphf<std::pair<uint64_t, uint64_t>, Custom_Pair_Hasher>(nelem, eqIdPair,
																				nthreads, gammaFactor);
	
	
	console -> info("BooPHF constructed perfect hash for {} keys; total memory = {} MB, bits/elem : {}\n",
					nelem, (mph -> totalBitSize() / 8) / (1024 * 1024), (double)(mph -> totalBitSize()) / nelem);


	auto t_end = time(nullptr);
	console -> info("Building the MPH table took time {} seconds.", t_end - t_start);
}





template <typename qf_obj, typename key_obj>
void ColoredDbg<qf_obj, key_obj>::
	build_color_classes_and_redirection_table(ColoredDbg<qf_obj, key_obj> &cdbg1, ColoredDbg<qf_obj, key_obj> &cdbg2)
{
	auto t_start = time(nullptr);

	console -> info("At color-class building (bitvectors concatenation) phase. Time = {}\n",
					time(nullptr) - start_time_);


	// Current buckets are required to be initialized with invalid values,
	// for the loading of the first BitVectorRRR file from each CdBG.
	uint64_t serialId = 0, currBucket1 = cdbg1.get_eq_class_file_count() + 1,
				currBucket2 = cdbg2.get_eq_class_file_count() + 1;
	
	BitVectorRRR bitVec1, bitVec2;


	for(auto idPair = eqIdPair.begin(); idPair != eqIdPair.end(); ++idPair)
	{
		uint64_t reqBucket1 = (idPair -> first ? (idPair -> first - 1) / mantis::NUM_BV_BUFFER : 0),
					reqBucket2 = (idPair -> second ? (idPair -> second - 1) / mantis::NUM_BV_BUFFER : 0);

		if(reqBucket1 != currBucket1)
		{
			sdsl::load_from_file(bitVec1, cdbg1.get_eq_class_files()[reqBucket1]);
			currBucket1 = reqBucket1;

			console -> info("Mantis 1: loaded one bitvectorRRR from file {}. Time = {}.",
							cdbg1.get_eq_class_files()[currBucket1], time(nullptr) - start_time_);
		}

		if(reqBucket2 != currBucket2)
		{
			sdsl::load_from_file(bitVec2, cdbg2.get_eq_class_files()[reqBucket2]);
			currBucket2 = reqBucket2;

			console -> info("Mantis 2: loaded one bitvectorRRR from file {}. Time = {}.",
							cdbg2.get_eq_class_files()[currBucket2], time(nullptr) - start_time_);
		}


		mphRedirect[mph -> lookup(*idPair)] = ++serialId;

		BitVector mergedColorClass(num_samples);
		concat(bitVec1, cdbg1.get_num_samples(), idPair -> first,
				bitVec2, cdbg2.get_num_samples(), idPair -> second, mergedColorClass);

		add_bitvector(mergedColorClass, serialId - 1);

		// Serialization and disk-write if required
		if(serialId % mantis::NUM_BV_BUFFER == 0)
		{
			console -> info("Serializing bitvector buffer with {} color-classes.", serialId);

			bv_buffer_serialize(serialId);
		}
	}


	// Serialize the bv buffer last time if needed
	if (serialId % mantis::NUM_BV_BUFFER > 0)
	{
		console -> info("Serializing bitvector buffer with {} color-classes.", serialId);
		bv_buffer_serialize(serialId);
	}


	auto t_end = time(nullptr);
	console -> info("Color-class building phase took time {} seconds.", t_end - t_start);
}





template <typename qf_obj, typename key_obj>
void ColoredDbg<qf_obj, key_obj>::
	build_cqf(ColoredDbg<qf_obj, key_obj> &cdbg1, ColoredDbg<qf_obj, key_obj> &cdbg2)
{
	auto t_start = time(nullptr);

	console -> info("At CQFs merging phase. Time = {}\n", time(nullptr) - start_time_);


	const CQF<key_obj> *cqf1 = cdbg1.get_cqf(), *cqf2 = cdbg2.get_cqf();
	typename CQF<key_obj>::Iterator it1 = cqf1 -> begin(), it2 = cqf2 -> begin(),
									walkBehindIterator1, walkBehindIterator2, walkBehindIteratorOut;
	uint64_t cqfPosition1 = 0, cqfPosition2 = 0, cqfOutPosition = 0;

	uint64_t kmerCount = 0, equalKmerCount = 0;
	uint64_t kmer1, kmer2, kmer, eqClass1, eqClass2;
	key_obj cqfEntry1, cqfEntry2;

	
	if(it1.done() || it2.done())
		console -> error("One or more CQF iterator(s) already at end position before starting walk.");



	
	if(!it1.done())
		cqfEntry1 = it1.get_cur_hash(), cqfPosition1++;
	
	if(!it2.done())
		cqfEntry2 = it2.get_cur_hash(), cqfPosition2++;


	while(!it1.done() || !it2.done())
	{
		if(!it1.done())
			kmer1 = cqfEntry1.key;
		
		if(!it2.done())
			kmer2 = cqfEntry2.key;


		if(it1.done())
		{
			kmer = kmer2;

			eqClass1 = 0, eqClass2 = cqfEntry2.count;
			advance_iterator_window(it2, cqfPosition2, cqfEntry2, walkBehindIterator2, cqf2);// ++it2;
		}
		else if(it2.done())
		{
			kmer = kmer1;
			
			eqClass1 = cqfEntry1.count, eqClass2 = 0;
			advance_iterator_window(it1, cqfPosition1, cqfEntry1, walkBehindIterator1, cqf1);// ++it1;
		}
		else if(kmer1 < kmer2)
		{
			kmer = kmer1;
			
			eqClass1 = cqfEntry1.count, eqClass2 = 0;
			advance_iterator_window(it1, cqfPosition1, cqfEntry1, walkBehindIterator1, cqf1);// ++it1;
		}
		else if(kmer2 < kmer1)
		{
			kmer = kmer2;

			eqClass1 = 0, eqClass2 = cqfEntry2.count;
			advance_iterator_window(it2, cqfPosition2, cqfEntry2, walkBehindIterator2, cqf2);// ++it2;
		}
		else
		{
			kmer = kmer1;

			eqClass1 = cqfEntry1.count, eqClass2 = cqfEntry2.count;
			advance_iterator_window(it1, cqfPosition1, cqfEntry1, walkBehindIterator1, cqf1),// ++it1;
			advance_iterator_window(it2, cqfPosition2, cqfEntry2, walkBehindIterator2, cqf2);// ++it2;

			equalKmerCount++; // for debugging purpose(s)
		}


		uint64_t colorId = mphRedirect[mph -> lookup(std::make_pair(eqClass1, eqClass2))];
		add_kmer(kmer, colorId, cqfOutPosition, walkBehindIteratorOut);
		// abundance[colorId]++;
		kmerCount++;

		if(kmerCount % PROGRESS_STEP == 0)
			console -> info("Kmers merged: {}M, time: {}", kmerCount * 10 / PROGRESS_STEP,
							time(nullptr) - start_time_);
	}


	console -> info("Total kmers merged: {}. Time: {}", kmerCount, time(nullptr) - start_time_);


	auto t_end = time(nullptr);
	console -> info("Merging the CQFs took time {} seconds.", t_end - t_start);
}











template <typename qf_obj, typename key_obj>
void ColoredDbg<qf_obj, key_obj>::
	serialize_cqf_and_abundance_dist(ColoredDbg<qf_obj, key_obj> &cdbg1, ColoredDbg<qf_obj, key_obj> &cdbg2)
{
	// Serialize the CQF
	if(dbg_alloc_flag == MANTIS_DBG_IN_MEMORY)
		dbg.serialize(prefix + mantis::CQF_FILE);
	else
		dbg.close();


	// Serialize the bv buffer last time if needed.
	// Done at the color-class building phase.


	// TODO: Move the list concatenation phase to constructor, with introduction of a new field;
	// just keep the disk-write part here.

	// Serialize the sample id map
	std::ofstream outputFile(prefix + mantis::SAMPLEID_FILE);

	for(auto sampleID: cdbg1.get_sample_id_map())
		outputFile << sampleID.first << " " << sampleID.second << "\n";

	for(auto sampleID: cdbg2.get_sample_id_map())
		outputFile << cdbg1.get_num_samples() + sampleID.first << " " << sampleID.second << "\n";

	outputFile.close();


	// Dump the abundance distribution of the equivalence / color classes.
	// if (flush_eqclass_dis)
	// {
	// 	const char OUTPUT_ABUNDANCE_DIST_FILE[] = "eqclass_dist.lst";
	// 	std::ofstream output(prefix + OUTPUT_ABUNDANCE_DIST_FILE);

	// 	for(uint64_t i = 1; i < abundance.size(); ++i)
	// 		output << i << " " << abundance[i] << "\n";
		
	// 	output.flush();
	// 	output.close();

	// 	abundance.clear();
	// 	abundance.shrink_to_fit();
	// }
}





template <typename qf_obj, typename key_obj>
void ColoredDbg<qf_obj, key_obj>::
	init_disk_buckets(ColoredDbg<qf_obj, key_obj> &cdbg1, ColoredDbg<qf_obj, key_obj> &cdbg2,
					std::vector<std::vector<std::ofstream>> &diskBucket)
{
	const uint64_t fileCount1 = cdbg1.get_eq_class_file_count(),
					fileCount2 = cdbg2.get_eq_class_file_count();

	console -> info("Initializing {} x {} disk-buckets.", fileCount1, fileCount2);

	diskBucket.resize(fileCount1 + 1),
	bucketSize.resize(fileCount1 + 1),
	cumulativeBucketSize.resize(fileCount1 + 1),
	MPH.resize(fileCount1 + 1);

	for(uint64_t i = 0; i <= fileCount1; ++i)
	{
		diskBucket[i].resize(fileCount2 + 1),
		bucketSize[i].resize(fileCount2 + 1),
		cumulativeBucketSize[i].resize(fileCount2 + 1),
		MPH[i].resize(fileCount2 + 1);
		
		for(uint64_t j = 0; j <= fileCount2; ++j)
			diskBucket[i][j] = std::ofstream(prefix + TEMP_DIR + EQ_ID_PAIRS_FILE + "_" +
											std::to_string(i) + "_" + std::to_string(j)),
			bucketSize[i][j] = 0,
			cumulativeBucketSize[i][j] = 0,
			MPH[i][j] = NULL;
	}

	console -> info("{} x {} disk-buckets initialized.", fileCount1, fileCount2);
}





template <typename qf_obj, typename key_obj>
void ColoredDbg<qf_obj, key_obj>::
	add_color_id_pair(const uint64_t colorID1, const uint64_t colorID2,
						std::vector<std::vector<std::ofstream>> &diskBucket)
{
	// TODO: Add faster file-write mechanism.

	const uint64_t row = (colorID1 ? (colorID1 - 1) / mantis::NUM_BV_BUFFER + 1 : 0),
					col = (colorID2 ? (colorID2 - 1) / mantis::NUM_BV_BUFFER + 1 : 0);

	diskBucket[row][col] << colorID1 << " " << colorID2 << "\n";
}





template <typename qf_obj, typename key_obj>
uint64_t ColoredDbg<qf_obj, key_obj>::
	fill_disk_buckets(ColoredDbg<qf_obj, key_obj> &cdbg1, ColoredDbg<qf_obj, key_obj> &cdbg2,
							std::unordered_set<std::pair<uint64_t, uint64_t>, Custom_Pair_Hasher> &sampledPairs)
{
	auto t_start = time(nullptr);

	console -> info("Writing all the color-id pairs to disk-files ({}). Time-stamp = {}",
					TEMP_DIR + EQ_ID_PAIRS_FILE + std::string("_X_Y"), time(nullptr) - start_time_);

	
	std::vector<std::vector<std::ofstream>> diskBucket;

	init_disk_buckets(cdbg1, cdbg2, diskBucket);


	uint64_t writtenPairsCount = 0;

	for(auto p = sampledPairs.begin(); p != sampledPairs.end(); ++p)
		add_color_id_pair(p -> first, p -> second, diskBucket);

	writtenPairsCount = sampledPairs.size();

	console -> info("Wrote the {} sampled pairs to disk. Time-stamp = {}", writtenPairsCount,
					time(nullptr) - start_time_);



	console -> info("Now iterating over the CQFs for the rest of the id-pairs.");

	const CQF<key_obj> *cqf1 = cdbg1.get_cqf(), *cqf2 = cdbg2.get_cqf();
	typename CQF<key_obj>::Iterator it1 = cqf1 -> begin(), it2 = cqf2 -> begin(),
									walkBehindIterator1, walkBehindIterator2;
	uint64_t cqfPosition1 = 0, cqfPosition2 = 0;

	uint64_t kmerCount = 0;
	uint64_t kmer1, kmer2;// eqClass1, eqClass2;
	std::pair<uint64_t, uint64_t> idPair;
	key_obj cqfEntry1, cqfEntry2;


	if(it1.done() || it2.done())
		console -> error("One or more CQF iterator(s) already at end position before starting walk.");


	if(!it1.done())
		cqfEntry1 = it1.get_cur_hash(), cqfPosition1++;
	
	if(!it2.done())
		cqfEntry2 = it2.get_cur_hash(), cqfPosition2++;


	while(!it1.done() || !it2.done())
	{
		if(!it1.done())
			kmer1 = cqfEntry1.key;
		
		if(!it2.done())
			kmer2 = cqfEntry2.key;


		// eqClassX = 0 implies absence in CdBG X.
		if(it1.done())
		{
			// eqClass1 = 0, eqClass2 = cqfEntry2.count;
			idPair.first = 0, idPair.second = cqfEntry2.count;
			advance_iterator_window(it2, cqfPosition2, cqfEntry2, walkBehindIterator2, cqf2);// ++it2;
		}
		else if(it2.done())
		{
			// eqClass1 = cqfEntry1.count, eqClass2 = 0;
			idPair.first = cqfEntry1.count, idPair.second = 0;
			advance_iterator_window(it1, cqfPosition1, cqfEntry1, walkBehindIterator1, cqf1);// ++it1;
		}
		else if(kmer1 < kmer2)
		{
			// eqClass1 = cqfEntry1.count, eqClass2 = 0;
			idPair.first = cqfEntry1.count, idPair.second = 0;
			advance_iterator_window(it1, cqfPosition1, cqfEntry1, walkBehindIterator1, cqf1);// ++it1;
		}
		else if(kmer2 < kmer1)
		{
			// eqClass1 = 0, eqClass2 = cqfEntry2.count;
			idPair.first = 0, idPair.second = cqfEntry2.count;
			advance_iterator_window(it2, cqfPosition2, cqfEntry2, walkBehindIterator2, cqf2);// ++it2;
		}
		else
		{
			// eqClass1 = cqfEntry1.count, eqClass2 = cqfEntry2.count;
			idPair.first = cqfEntry1.count, idPair.second = cqfEntry2.count;
			advance_iterator_window(it1, cqfPosition1, cqfEntry1, walkBehindIterator1, cqf1),// ++it1;
			advance_iterator_window(it2, cqfPosition2, cqfEntry2, walkBehindIterator2, cqf2);// ++it2;
		}


		kmerCount++;

		if(sampledPairs.find(idPair) == sampledPairs.end())
		{
			// output << eqClass1 << " " << eqClass2 << "\n";
			add_color_id_pair(idPair.first, idPair.second, diskBucket);
			writtenPairsCount++;
		}


		if(kmerCount % PROGRESS_STEP == 0)
			console -> info("Observed count of distinct k-mers: {}M, written id-pairs to disk: {}. Time-stamp = {}",
							kmerCount * 10 / PROGRESS_STEP, writtenPairsCount, time(nullptr) - start_time_);
	}


	console -> info("Distinct kmers found {}, id pairs written to disk {}. Time-stamp = {}",
					kmerCount, writtenPairsCount, time(nullptr) - start_time_);

	// output.flush();
	// output.close();

	// Flush and close the disk-bucket files.
	for(int i = 0; i <= cdbg1.get_eq_class_file_count(); ++i)
		for(int j = 0; j <= cdbg1.get_eq_class_file_count(); ++j)
		{
			diskBucket[i][j].flush();
			diskBucket[i][j].close();
		}

	
	auto t_end = time(nullptr);
	console -> info("Filling up the disk-buckets with color-id pairs took time {} seconds.", t_end - t_start);


	return kmerCount;
}





template <typename qf_obj, typename key_obj>
uint64_t ColoredDbg<qf_obj, key_obj>::
	filter_disk_buckets(ColoredDbg<qf_obj, key_obj> &cdbg1, ColoredDbg<qf_obj, key_obj> &cdbg2)
{
	if(!system(NULL))
	{
		console -> error("Command processor does not exist.");
		exit(1);
	}


	auto t_start = time(nullptr);

	console -> info("Filtering out the unique eq-id pairs from files {} with {} threads. Time-stamp = {}",
					TEMP_DIR + EQ_ID_PAIRS_FILE + "_X_Y", threadCount, time(nullptr) - start_time_);


	uint64_t colorClassCount = 0;
	const uint64_t fileCount1 = cdbg1.get_eq_class_file_count(),
					fileCount2 = cdbg2.get_eq_class_file_count();

	
	for(int i = 0; i <= fileCount1; ++i)
		for(int j = 0; j <= fileCount2; ++j)
		{
			std::string diskBucket = prefix + TEMP_DIR + EQ_ID_PAIRS_FILE +
									"_" + std::to_string(i) + "_" + std::to_string(j);

			std::string sysCommand = "sort -u";
			// TODO: Consult professor on parallelization and memory-usage details.
			sysCommand += " --parallel=" + std::to_string(threadCount);
			sysCommand += " -S " + std::to_string(maxMemoryForSort) + "G";
			sysCommand += " -o " + diskBucket + " " + diskBucket;

			console -> info("System command used:\n{}", sysCommand);

			system(sysCommand.c_str());

			
			std::string lineCountFile = prefix + TEMP_DIR + ID_PAIR_COUNT_FILE;
			sysCommand = "wc -l " + diskBucket + " | egrep -o \"[0-9]+ \" > " + lineCountFile;
			system(sysCommand.c_str());

			std::ifstream lineCount(lineCountFile);

			lineCount >> bucketSize[i][j];
			lineCount.close();


			colorClassCount += bucketSize[i][j];

			console -> info("Filtered {} unique color-id pairs of disk-bucket {}. Time-stamp = {}",
							bucketSize[i][j], diskBucket, time(nullptr) - start_time_);
		}


	console -> info("Count of unique color-id pairs = {}. Time-stamp = {}",
					colorClassCount, time(nullptr) - start_time_);

	

	auto t_end = time(nullptr);
	console -> info("Filtering the unique color-id pairs took time {} seconds.", t_end - t_start);

	return colorClassCount;
}





template <typename qf_obj, typename key_obj>
void ColoredDbg<qf_obj, key_obj>::
	build_MPH_tables(ColoredDbg<qf_obj, key_obj> &cdbg1, ColoredDbg<qf_obj, key_obj> &cdbg2)
{	
	auto t_start = time(nullptr);

	console -> info("Building an MPH (Minimal Perfect Hash) table per disk-bucket. Time-stamp = {}",
					time(nullptr) - start_time_);


	const uint64_t fileCount1 = cdbg1.get_eq_class_file_count(),
					fileCount2 = cdbg2.get_eq_class_file_count();

	uint64_t mphMemory = 0;

	
	for(uint64_t i = 0; i <= fileCount1; ++i)
		for(uint64_t j = 0; j <= fileCount2; ++j)
			if(!bucketSize[i][j])
				console -> info("Bucket ({}, {}) is empty. No MPH table is built.", i, j);
			else
			{
				std::vector<std::pair<uint64_t, uint64_t>> colorIdPair;
				colorIdPair.reserve(bucketSize[i][j]);


				console -> info("Loading the unique color-id pairs from bucket ({}, {}) into memory. Time-stamp = {}",
								i, j, time(nullptr) - start_time_);

				std::ifstream input(prefix + TEMP_DIR + EQ_ID_PAIRS_FILE +
									"_" + std::to_string(i) + "_" + std::to_string(j));
				uint64_t id1, id2;

				while(input >> id1 >> id2)
					colorIdPair.emplace_back(id1, id2);

				input.close();

				console -> info("Loaded {} color-id pairs into memory. Time-stamp = {}", colorIdPair.size(),
								time(nullptr) - start_time_);

				
				console -> info("Constructing a BooPHF with {} elements using {} threads. Time-stamp = {}",
								colorIdPair.size(), threadCount, time(nullptr) - start_time_);

				// lowest bit/elem is achieved with gamma=1, higher values lead to larger mphf but faster construction/query
				double gammaFactor = 2.0;	// gamma = 2 is a good tradeoff (leads to approx 3.7 bits/key )

				// build the mphf
				MPH[i][j] = new boomphf::mphf<std::pair<uint64_t, uint64_t>, Custom_Pair_Hasher>(colorIdPair.size(),
																								colorIdPair,
																								threadCount,
																								gammaFactor);
				mphMemory += MPH[i][j] -> totalBitSize();

				
				console -> info("For bucket ({}, {}) BooPHF constructed perfect hash for {} keys; total memory = {} MB, bits/elem : {}. Time-stamp = {}\n",
								i, j, colorIdPair.size(), (MPH[i][j] -> totalBitSize() / 8) / (1024 * 1024),
								(double)(MPH[i][j] -> totalBitSize()) / colorIdPair.size(), time(nullptr) - start_time_);

				
				colorIdPair.clear();
				colorIdPair.shrink_to_fit();
			}

	console -> info("Total memory consumed by all the MPH tables = {} MB.", (mphMemory / 8) / (1024 * 1024));
	

	auto t_end = time(nullptr);
	console -> info("Building the MPH tables took time {} seconds.", t_end - t_start);
}





template <typename qf_obj, typename key_obj>
void ColoredDbg<qf_obj, key_obj>::
	build_color_class_table(ColoredDbg<qf_obj, key_obj> &cdbg1, ColoredDbg<qf_obj, key_obj> &cdbg2)
{
	auto t_start = time(nullptr);

	console -> info("At color-class building (bitvectors concatenation) phase. Time = {}\n",
					time(nullptr) - start_time_);


					
	const uint64_t fileCount1 = cdbg1.get_eq_class_file_count(), fileCount2 = cdbg2.get_eq_class_file_count();

	uint64_t writtenPairsCount = 0;
	BitVectorRRR bitVec1, bitVec2;


	for(uint64_t i = 0; i <= fileCount1; ++i)
	{
		if(i > 0)
		{
			sdsl::load_from_file(bitVec1, cdbg1.get_eq_class_files()[i - 1]);
			console -> info("Mantis 1: loaded one bitvectorRRR from file {}. Time = {}.",
							cdbg1.get_eq_class_files()[i - 1], time(nullptr) - start_time_);
		}

		for(uint64_t j = 0; j <= fileCount2; ++j)
		{
			if(j > 0)
			{
				sdsl::load_from_file(bitVec2, cdbg2.get_eq_class_files()[j - 1]);
				console -> info("Mantis 2: loaded one bitvectorRRR from file {}. Time = {}.",
								cdbg2.get_eq_class_files()[j - 1], time(nullptr) - start_time_);
			}

			console -> info("At bucket ({}, {}), size = {}. Time-stamp = {}", i, j, bucketSize[i][j],
							time(nullptr) - start_time_);


			cumulativeBucketSize[i][j] = writtenPairsCount;


			if(bucketSize[i][j])
			{
				uint64_t queueCount = ((cumulativeBucketSize[i][j] + bucketSize[i][j] - 1) / mantis::NUM_BV_BUFFER)
										- (cumulativeBucketSize[i][j] / mantis::NUM_BV_BUFFER) + 1;

				// std::vector<std::deque<std::pair<uint64_t, uint64_t>>> writeQueue;
				std::vector<std::vector<std::pair<uint64_t, uint64_t>>> writeQueue;

				writeQueue.resize(queueCount);

				uint64_t remPairsCount = bucketSize[i][j];
				for(uint64_t k = 0; k < queueCount; ++k)
				{
					uint64_t queueLen = std::min(remPairsCount, mantis::NUM_BV_BUFFER -
							(cumulativeBucketSize[i][j] + (bucketSize[i][j] - remPairsCount)) % mantis::NUM_BV_BUFFER);
					writeQueue[k].reserve(queueLen);

					remPairsCount -= queueLen;
				}



				std::ifstream input(prefix + TEMP_DIR + EQ_ID_PAIRS_FILE +
									"_" + std::to_string(i) + "_" + std::to_string(j));
				std::pair<uint64_t, uint64_t> idPair;

				while(input >> idPair.first >> idPair.second)
				{
					uint64_t queueIdx = ((cumulativeBucketSize[i][j] + MPH[i][j] -> lookup(idPair)) / mantis::NUM_BV_BUFFER)
										- (cumulativeBucketSize[i][j] / mantis::NUM_BV_BUFFER);
					writeQueue[queueIdx].push_back(idPair);
				}

				input.close();



				for(uint64_t k = 0; k < queueCount; ++k)
				{
					for(auto it = writeQueue[k].begin(); it != writeQueue[k].end(); ++it)
					{
						uint64_t colorID = cumulativeBucketSize[i][j] + MPH[i][j] -> lookup(*it) + 1;

						BitVector colorClass(num_samples);
						concat(bitVec1, cdbg1.get_num_samples(), it -> first,
								bitVec2, cdbg2.get_num_samples(), it -> second, colorClass);

						add_bitvector(colorClass, colorID - 1), writtenPairsCount++;

						if(writtenPairsCount % mantis::NUM_BV_BUFFER == 0)
						{
							console -> info("Serializing bitvector buffer with {} color-classes.", writtenPairsCount);

							bv_buffer_serialize(writtenPairsCount);
						}
					}

					
					writeQueue[k].clear();
					writeQueue[k].shrink_to_fit();
				}


				writeQueue.clear();
				writeQueue.shrink_to_fit();
			}
		}
	}


	// Serialize the bv buffer last time if needed
	if (writtenPairsCount % mantis::NUM_BV_BUFFER > 0)
	{
		console -> info("Serializing bitvector buffer with {} color-classes.", writtenPairsCount);
		bv_buffer_serialize(writtenPairsCount);
	}


	auto t_end = time(nullptr);
	console -> info("Color-class building phase took time {} seconds.", t_end - t_start);
}





template <typename qf_obj, typename key_obj>
uint64_t ColoredDbg<qf_obj, key_obj>:: get_color_id(const std::pair<uint64_t, uint64_t> &idPair)
{
	const uint64_t row = (idPair.first ? (idPair.first - 1) / mantis::NUM_BV_BUFFER + 1 : 0),
					col = (idPair.second ? (idPair.second - 1) / mantis::NUM_BV_BUFFER + 1 : 0);


	return cumulativeBucketSize[row][col] + MPH[row][col] -> lookup(idPair) + 1;
}





template <typename qf_obj, typename key_obj>
void ColoredDbg<qf_obj, key_obj>::
	build_CQF(ColoredDbg<qf_obj, key_obj> &cdbg1, ColoredDbg<qf_obj, key_obj> &cdbg2)
{
	auto t_start = time(nullptr);

	console -> info("At CQFs merging phase. Time = {}\n", time(nullptr) - start_time_);


	const CQF<key_obj> *cqf1 = cdbg1.get_cqf(), *cqf2 = cdbg2.get_cqf();
	typename CQF<key_obj>::Iterator it1 = cqf1 -> begin(), it2 = cqf2 -> begin(),
									walkBehindIterator1, walkBehindIterator2, walkBehindIteratorOut;
	uint64_t cqfPosition1 = 0, cqfPosition2 = 0, cqfOutPosition = 0;

	uint64_t kmerCount = 0, equalKmerCount = 0;
	uint64_t kmer1, kmer2, kmer;//, eqClass1, eqClass2;
	std::pair<uint64_t, uint64_t> idPair;
	key_obj cqfEntry1, cqfEntry2;

	
	if(it1.done() || it2.done())
		console -> error("One or more CQF iterator(s) already at end position before starting walk.");



	
	if(!it1.done())
		cqfEntry1 = it1.get_cur_hash(), cqfPosition1++;
	
	if(!it2.done())
		cqfEntry2 = it2.get_cur_hash(), cqfPosition2++;


	while(!it1.done() || !it2.done())
	{
		if(!it1.done())
			kmer1 = cqfEntry1.key;
		
		if(!it2.done())
			kmer2 = cqfEntry2.key;


		if(it1.done())
		{
			kmer = kmer2;

			idPair.first = 0, idPair.second = cqfEntry2.count;
			advance_iterator_window(it2, cqfPosition2, cqfEntry2, walkBehindIterator2, cqf2);// ++it2;
		}
		else if(it2.done())
		{
			kmer = kmer1;
			
			idPair.first = cqfEntry1.count, idPair.second = 0;
			advance_iterator_window(it1, cqfPosition1, cqfEntry1, walkBehindIterator1, cqf1);// ++it1;
		}
		else if(kmer1 < kmer2)
		{
			kmer = kmer1;
			
			idPair.first = cqfEntry1.count, idPair.second = 0;
			advance_iterator_window(it1, cqfPosition1, cqfEntry1, walkBehindIterator1, cqf1);// ++it1;
		}
		else if(kmer2 < kmer1)
		{
			kmer = kmer2;

			idPair.first = 0, idPair.second = cqfEntry2.count;
			advance_iterator_window(it2, cqfPosition2, cqfEntry2, walkBehindIterator2, cqf2);// ++it2;
		}
		else
		{
			kmer = kmer1;

			idPair.first = cqfEntry1.count, idPair.second = cqfEntry2.count;
			advance_iterator_window(it1, cqfPosition1, cqfEntry1, walkBehindIterator1, cqf1),// ++it1;
			advance_iterator_window(it2, cqfPosition2, cqfEntry2, walkBehindIterator2, cqf2);// ++it2;

			equalKmerCount++; // for debugging purpose(s)
		}


		// uint64_t colorId = mphRedirect[mph -> lookup(std::make_pair(eqClass1, eqClass2))]; 
		add_kmer(kmer, get_color_id(idPair), cqfOutPosition, walkBehindIteratorOut);
		// abundance[colorId]++;
		kmerCount++;

		if(kmerCount % PROGRESS_STEP == 0)
			console -> info("Kmers merged: {}M, time: {}", kmerCount * 10 / PROGRESS_STEP,
							time(nullptr) - start_time_);
	}


	console -> info("Total kmers merged: {}. Time: {}", kmerCount, time(nullptr) - start_time_);


	auto t_end = time(nullptr);
	console -> info("Merging the CQFs took time {} seconds.", t_end - t_start);
}





template <typename qf_obj, typename key_obj>
uint64_t ColoredDbg<qf_obj, key_obj>::
	merge_2(ColoredDbg<qf_obj, key_obj> &cdbg1, ColoredDbg<qf_obj, key_obj> &cdbg2)
{
	auto t_start = time(nullptr);
	console -> info ("Merge starting. Time = {}\n", time(nullptr) - start_time_);


	// Make the temporary directory if it doesn't exist.
	std::string tempDir = prefix + TEMP_DIR;

	if(!mantis::fs::DirExists(tempDir.c_str()))
		mantis::fs::MakeDir(tempDir.c_str());
	
	// Check to see if the temporary directory exists now.
	if(!mantis::fs::DirExists(tempDir.c_str()))
	{
		console->error("Temporary directory {} could not be successfully created.", tempDir);
		exit(1);
	}


	const uint64_t samplePairCount = 1000000;
	std::unordered_set<std::pair<uint64_t, uint64_t>, Custom_Pair_Hasher> sampledPairs;
	
	sample_eq_id_pairs(cdbg1, cdbg2, mantis::SAMPLE_SIZE, samplePairCount, sampledPairs);
	uint64_t kmerCount = fill_disk_buckets(cdbg1, cdbg2, sampledPairs);
	
	sampledPairs.clear();

	uint64_t colorClassCount = filter_disk_buckets(cdbg1, cdbg2);

	build_MPH_tables(cdbg1, cdbg2);

	bv_buffer = BitVector(mantis::NUM_BV_BUFFER * num_samples);
	build_color_class_table(cdbg1, cdbg2);
	bv_buffer = BitVector(0);

	initialize_CQF(cdbg1.get_cqf() -> keybits(), cdbg1.get_cqf() -> hash_mode(), cdbg1.get_cqf() -> seed(), 
					kmerCount);
	build_CQF(cdbg1, cdbg2);


	// eqIdPair.reserve(colorClassCount);

	// load_unique_id_pairs();
	// block_sort();
	// build_mph_table();

	// mphRedirect.resize(colorClassCount);
	// bv_buffer = BitVector(mantis::NUM_BV_BUFFER * num_samples);

	// build_color_classes_and_redirection_table(cdbg1, cdbg2);

	// eqIdPair.clear();
	// eqIdPair.shrink_to_fit();
	// bv_buffer = BitVector(0);

	// // abundance.resize(colorClassCount + 1, 0);

	// initialize_CQF(cdbg1.get_cqf() -> keybits(), cdbg1.get_cqf() -> hash_mode(), cdbg1.get_cqf() -> seed(), 
	// 				kmerCount);
	// build_cqf(cdbg1, cdbg2);

	// delete mph;
	// mphRedirect.clear();
	// mphRedirect.shrink_to_fit();


	// Remove the temporary directory
	// std::string sysCommand = "rm -rf " + tempDir;
	// system(sysCommand.c_str());


	auto t_end = time(nullptr);
	console -> info("Merge completed. Total time taken is {} seconds.", t_end - t_start);

	return colorClassCount;
}

#endif
