/*
 * ============================================================================
 *
 *         Author:  Jamshed Khan, jamshed@cs.umd.edu
 *   Organization:  University of Maryland
 *
 * ============================================================================
 */

#ifndef _CDBG_MERGER_H_
#define _CDBG_MERGER_H_

#include "coloreddbg.h"
#include "BooPHF.h"

#include "mstMerger.h"

template <class qf_obj, class key_obj>
class CdBG_Merger
{
	public:
        CdBG_Merger(ColoredDbg<qf_obj, key_obj> &cdbg1, ColoredDbg<qf_obj, key_obj> &cdbg2,
                    ColoredDbg<qf_obj, key_obj> &cdbgOut);

        void set_console(spdlog::logger* c) { console = c; }

		// Sets the number of processor-threads to be used at the intermediate steps of 
		// unique id-pairs filtering and MPH building.
		inline void set_thread_count(uint threadNum) { threadCount = threadNum; }

		// Merges two Colored dBG 'cdbg1' and 'cdbg2' into the colored dBG 'cdbg'
		// (all the CdBG's are private members).
		void merge();


	private:
		// k-mer count in progress display.
		const static uint64_t PROGRESS_STEP = 10000000;

		// CQF-window size to keep in memory.
		const static uint64_t ITERATOR_WINDOW_SIZE = 4096;

		// Count of popular color-id pairs to be sampled
		const static uint64_t SAMPLE_PAIR_COUNT = std::min((uint64_t)1000000, mantis::NUM_BV_BUFFER);

		// Name of the temporary working directory at disk; will be present
		// temporarily inside the output directory.
		const std::string TEMP_DIR = std::string("temp/");

		// Name of the temporary list of color-id pairs.
		const std::string EQ_ID_PAIRS_FILE = std::string("color-id-pairs");

		// Name of the temporary file to contain count of distinct color-id pairs.
		const std::string ID_PAIR_COUNT_FILE = std::string("color-id-pairs-count");

        // First input colored dBG.
        ColoredDbg<qf_obj, key_obj> cdbg1;

        // Second input colored dBG.
        ColoredDbg<qf_obj, key_obj> cdbg2;

        // Output (merged) colored dBG.
        ColoredDbg<qf_obj, key_obj> cdbg;

        // Console to display messages.
        spdlog::logger* console;

        // Initial time-stamp of object creation.
        std::time_t start_time_;
		
		// Number of processor-threads to be used at the intermediate steps of  unique
		// color-id pairs filtering and MPH's building.
		uint threadCount = 1;

        // Utility information to display at the result summary.
        uint64_t colorCount1 = 0, colorCount2 = 0;

        // Blocks for minimizers in the output CDBG
		std::vector<uint64_t> minimizerBlocks;
		std::vector<std::unordered_map<uint64_t, std::pair<colorIdType , colorIdType>>> minimizerMap;
		uint64_t kbits;


	// Required to hash colo-id pair objects. Resorted to boost::hash_combine
		// instead of plain XOR hashing. For more explanation, consult
		// https://stackoverflow.com/questions/35985960/c-why-is-boosthash-combine-the-best-way-to-combine-hash-values
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

		// Hash-map type for color-id pairs.
		typedef std::unordered_map<std::pair<uint64_t, uint64_t>, uint64_t, Custom_Pair_Hasher> idPairMap_t;

		// Hash-map for the sampled (on abundance) color-id pairs.
		// Used as the form (pair -> abundance) earlier, and finally as (pair -> colorId).
		idPairMap_t sampledPairs;

		// Disk-bucket filestreams.
		std::vector<std::vector<std::ofstream>> diskBucket;

		// Stores the size (number of color-id pairs) at each disk-bucket.
		std::vector<std::vector<uint64_t>> bucketSize;

		// The (i, j)'th entry contains the total count of color-id pairs upto
		// bucket (i, j), exclusive, in row-major order.
		std::vector<std::vector<uint64_t>> cumulativeBucketSize;

		// MPH (Minimal Perfect Hash) function tables for each of the disk-buckets.
		std::vector<std::vector<boophf_t *>> MPH;



        // Advances the CQF iterator 'it', with keeping track of the 'step' count; and
		// fetches the next CQF-entry into 'cqfEntry' if the iterator 'it' is advanced
		// into a non-end position. Also, advances or initializes the iterator
		// 'walkBehindIterator' that trails 'it' by ITERATOR_WINDOW_SIZE.
		static void advance_iterator_window(typename CQF<key_obj>::Iterator &it, uint64_t &step, key_obj &cqfEntry,
											typename CQF<key_obj>::Iterator &walkBehindIterator,
											const CQF<key_obj> *cqf);

        // Samples 'SAMPLE_PAIR_COUNT' number of most abundant color-id pairs from the
		// first 'sampleKmerCount' distinct k-mers of the CdBGs 'cdbg1' and 'cdbg2',
		// into the map 'sampledPairs', which is of the format (pair -> abundance).
		uint64_t sample_color_id_pairs(uint64_t sampleKmerCount);

		// Initializes the disk-buckets, i.e. initializes the disk-files, MPH tables,
		// bucket sizes, cumulative size counts etc.
		inline void init_disk_buckets();

        // Gathers all the color-id pairs for all the distinct k-mers of the CdBGs
		// 'cdbg1' and 'cdbg2' into disk-files (or referred to as buckets hereafter),
		// where an id-pair goes to bucket_(i, j) if it reads color-class bitvectors
		// from bitvector_file_(i-1) of cdbg1 and bitvector_file_(j-1) of cdbg2; with
		// avoiding writes of sampled pairs from the map 'sampledPairs', sampled on
		// abundance. A bucket of the form (0, X) with X > 0 implies that, the id-pairs
		// present at this bucket are only of k-mers that are absent at cdbg1, present
		// at cdbg2, and read from the bitvector_file_(X - 1) of cdbg2. Buckets of the
		// form (X, 0) imply vice versa.
		// Returns the number of distinct k-mers present at the CdBGs cdg1 and cdbg2.
		uint64_t fill_disk_buckets(uint64_t startingBlock = 0);

		// Adds the color-class ID pair (colorID1, colorID2) to the appropriate disk
		// bucket; i.e. writes the pair into the file diskBucket[i][j] iff colorID1
		// has its color-class (bitvector) at bitvector_file_(i - 1) and colorID2 has
		// its color-class at bitvector_file_(j - 1).
		inline void add_color_id_pair(const uint64_t colorID1, const uint64_t colorID2,
										std::vector<std::vector<std::ofstream>> &diskBucket);

        // Filters the disk-buckets to contain only unique color-id pairs.
		// Returns the count of unique color-id pairs.
		uint64_t filter_disk_buckets();

		// Returns the maximum amount of memory (in GB) to use at the color-id pairs
		// filtering phase, specifically, in the 'sort -u' system command.
		uint64_t get_max_sort_memory();

		// Assign color-ids to the abudant color-id pairs present at the hash-map
		// 'sampledPairs' (of the form (pair -> abundance)), based on their abundance.
		// The hash-map 'sampledPairs' will contain the assigned color-ids for the
		// pairs after execution of this method, losing the abundance information.
		void assign_abundant_color_ids();

		// Builds an MPH (Minimal Perfect Hash) table for each disk-bucket.
		void build_MPH_tables();

        // Builds the output color-class bitvectors for the color-id pairs
		// that are not sampled on abundance, i,e. those that are written to disk.
		void build_color_class_table();

        // Builds the output color-class bitvectors for the abundant (sampled)
		// color-id pairs.
		void build_abundant_color_classes();

		// Concatenates the color-classes (BitVector objects) corresponding to the
		// color-IDs 'colorID1' and 'colorID2' respectively, from two CdBGs of sample
		// sizes 'colCount1' and 'colCount2' each. The partial BitVectorRRR objects
		// containing the BitVectors for the two color-IDs are 'bv1' and 'bv2'
		// respectively. The concatenated result BitVector is stored in 'resultVec'.
		void concat(const BitVectorRRR &bv1, const uint64_t colCount1, const uint64_t colorID1,
					const BitVectorRRR &bv2, const uint64_t colCount2, const uint64_t colorID2,
					BitVector &resultVec);

		// Serialize the bitvector buffer to disk, when the current constructed class
		// count is 'colorClsCount'.
		void bv_buffer_serialize(uint64_t colorClsCount);

		// Initializes the CQF that will contain 'kmerCount' number of k-mers after
		// mantii merge.
		void initialize_CQF(uint32_t keybits, qf_hashmode hashMode, uint32_t seed, uint64_t kmerCount);

		// Builds the output merged CQF.
		void build_CQF();

        // Given a color-id pair 'idPair', returns the newly assigned color-id of this
		// pair at the merged CdBG.
		inline uint64_t get_color_id(const std::pair<uint64_t, uint64_t> &idPair);

		// Add a k-mer 'kmer' with color-class ID 'colorID' into the CQF of this CdBG.
		void add_kmer(uint64_t kmer, uint64_t colorID, uint64_t &step,
						typename CQF<key_obj>::Iterator &walkBehindIterator);

		// Serializes the output CQF and sample-id mapping.
		void serialize_cqf_and_sampleid_list();

	// Builds the output color-class bitvectors for the color-id pairs.
	void store_color_pairs(ColoredDbg<qf_obj, key_obj> &cdbg1, ColoredDbg<qf_obj, key_obj> &cdbg2, uint64_t& numColorBuffers);
	uint64_t store_abundant_color_pairs(std::ofstream& output);
	void calc_mst_stats(ColoredDbg<qf_obj, key_obj> &cdbg1, ColoredDbg<qf_obj, key_obj> &cdbg2, std::string& dir1, std::string& dir2);

	void divide_minimizers_into_blocks(std::vector<uint64_t> &min1, std::vector<uint64_t> &min2);
};



template<typename qf_obj, typename key_obj>
CdBG_Merger<qf_obj, key_obj>::
    CdBG_Merger(ColoredDbg<qf_obj, key_obj> &cdbg1, ColoredDbg<qf_obj, key_obj> &cdbg2,
                ColoredDbg<qf_obj, key_obj> &cdbgOut)
{
    start_time_ = std::time(nullptr);

	if(cdbg1.get_eq_class_file_count() >= cdbg2.get_eq_class_file_count())
		this -> cdbg1 = cdbg1, this -> cdbg2 = cdbg2;
	else
	{
		this -> cdbg1 = cdbg2, this -> cdbg2 = cdbg1;

//		console -> info("Mantis indices are swapped.");
	}
	
    cdbg = cdbgOut;
}


template <typename qf_obj, typename key_obj>
inline void CdBG_Merger <qf_obj, key_obj> ::
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
void CdBG_Merger<qf_obj, key_obj>::divide_minimizers_into_blocks(
		std::vector<uint64_t> &min1,
		std::vector<uint64_t> &min2) {

	minimizerBlocks.resize(min1.size(), 0);
	uint64_t blockCnt{min1[0]+min2[0]}, block{0};
	for (auto i = 1; i < minimizerBlocks.size(); i++) {
		minimizerBlocks[i - 1] = block;
		if ((blockCnt + min1[i]+min2[i]) > block_kmer_threshold) {
			block++;
			blockCnt = 0;
		}
		blockCnt += (min1[i]+min2[i]);
	}
	minimizerBlocks[minimizerBlocks.size() - 1] = block;
	console->info("Initial # of blocks : {}", block+1);
	minimizerMap.resize(block+1);
}


template <typename qf_obj, typename key_obj>
uint64_t CdBG_Merger<qf_obj, key_obj>::
	sample_color_id_pairs(uint64_t sampleKmerCount)
{
	auto t_start = time(nullptr);

	uint64_t sampleCount = SAMPLE_PAIR_COUNT;
	console -> info("Sampling {} most-abundant color-id pairs from the first {} kmers. Time-stamp = {}.",
					sampleCount, sampleKmerCount, time(nullptr) - start_time_);
	

	uint64_t kmerCount = 0;
	idPairMap_t pairCount;

	kbits = cdbg1.get_current_cqf()->keybits();

	auto kmerMask = (1ULL << kbits) - 1;

	uint64_t curBlock{0};
	console->info("cdbg1.numBlocks={}, cdbg2.numBlocks={}", cdbg1.get_numBlocks(), cdbg2.get_numBlocks());
	uint64_t maxMinimizer{0}, minMinimizer{invalid};
	while(kmerCount < sampleKmerCount and
            (curBlock < cdbg1.get_numBlocks() or curBlock < cdbg2.get_numBlocks())) {
	    std::cerr << "Current Block="<<curBlock << "\n";
		uint64_t maxMinimizer1{0}, minMinimizer1{invalid},
		maxMinimizer2{0}, minMinimizer2{invalid};

		cdbg1.replaceCQFInMemory(curBlock);
		if (curBlock < cdbg1.get_numBlocks()) {
            const CQF<key_obj> *cqf1 = cdbg1.get_current_cqf();
            typename CQF<key_obj>::Iterator it1 = cqf1->begin();
            uint64_t cntr{0};
            uint64_t count{cqf1->dist_elts()};
			auto minmax = cdbg1.getMinMaxMinimizer(curBlock);
			minMinimizer1 = minmax.first;
			maxMinimizer1 = minmax.second;
			while (!it1.done()) {
				auto keyval = it1.get_cur_hash();
				auto key = hash_64i(keyval.key, kmerMask);
				auto minimizerPair = cdbg1.findMinimizer(key, kbits);
				std::vector<uint64_t> pairs{minimizerPair.first, minimizerPair.second};
				for (auto minimizer : pairs) {
					if (minimizer != invalid and minimizer >= minMinimizer1 and minimizer <= maxMinimizer1) {
						auto &map = minimizerMap[minimizer];
						if (map.find(keyval.key) == map.end()) {
							map[keyval.key] = std::make_pair(keyval.count, 0);
						} else if (map[keyval.key].first and map[keyval.key].first != keyval.count) {
							console->error("found a key with different color class in the second cqf.");
							std::exit(3);
						} else if (not map[keyval.key].first) {
							map[keyval.key].first = keyval.count;
						}
					}
				}
                ++it1;
                cntr++;
                if (cntr % 10000000 == 0) {
                    std::cerr << "\rFirst mantis: " << cntr << " out of " << count;
                }
            }
        }

		cdbg2.replaceCQFInMemory(curBlock);
		if (curBlock < cdbg2.get_numBlocks()) {
            const CQF<key_obj> *cqf2 = cdbg2.get_current_cqf();
            typename CQF<key_obj>::Iterator it2 = cqf2->begin();
            uint64_t cntr{0};
            uint64_t count{cqf2->dist_elts()};
			auto minmax = cdbg2.getMinMaxMinimizer(curBlock);
			minMinimizer2 = minmax.first;
			maxMinimizer2 = minmax.second;
			while (!it2.done()) {
				auto keyval = it2.get_cur_hash();
				auto key = hash_64i(keyval.key, kmerMask);
				auto minimizerPair = cdbg2.findMinimizer(key, kbits);
				std::vector<uint64_t> pairs{minimizerPair.first, minimizerPair.second};
				for (auto minimizer : pairs) {
					if (minimizer != invalid and minimizer >= minMinimizer2 and minimizer <= maxMinimizer2) {
						auto &map = minimizerMap[minimizer];
						if (map.find(keyval.key) == map.end()) {
							map[keyval.key] = std::make_pair(0, keyval.count);
						} else if (map[keyval.key].second and map[keyval.key].second != keyval.count) {
							console->error("found a key with different color class in the second cqf.");
							std::exit(3);
						} else if (not map[keyval.key].second) {
							map[keyval.key].second = keyval.count;
						}
					}
				}
                ++it2;
                cntr++;
                if (cntr % 10000000 == 0) {
                    std::cerr << "\rSecond mantis: " << cntr << " out of " << count;
                }
            }
        }

		minMinimizer = curBlock >= cdbg1.get_numBlocks() ? minMinimizer2 :
					   curBlock >= cdbg2.get_numBlocks() ? minMinimizer1 : std::min(minMinimizer1, minMinimizer2);
		maxMinimizer = curBlock >= cdbg1.get_numBlocks() ? maxMinimizer2 :
					   curBlock >= cdbg2.get_numBlocks() ? maxMinimizer1 : std::min(maxMinimizer1, maxMinimizer2);
		std::cerr << "\nMin minimizer=" << minMinimizer << " Max minimizer=" << maxMinimizer << "\n";
		curBlock++;
		uint64_t b = minMinimizer;
		while (b <= maxMinimizer) {
			if (not minimizerMap[b].empty()) {
				std::cerr << "\rminimizer " << b << " size=" << minimizerMap[b].size() << "     ";
				for (auto &kv : minimizerMap[b]) {
					pairCount[kv.second]++;
					kmerCount++;
				}
			}
			minimizerMap[b].clear();
			b++;
		}
		std::cerr << "\r";
		minMinimizer = maxMinimizer + 1;
	}
	console -> info("Sampled {} k-mers, color-classes found: {}. Time-stamp = {}.",
					kmerCount, pairCount.size(), time(nullptr) - start_time_);

	typedef std::pair<uint64_t, std::pair<uint64_t, uint64_t>> CountAndIdPair;
	std::priority_queue<CountAndIdPair, std::vector<CountAndIdPair>, std::greater<CountAndIdPair>> minPQ;

	for(auto p = pairCount.begin(); p != pairCount.end(); ++p)
		if(minPQ.size() < SAMPLE_PAIR_COUNT)
			minPQ.push(std::make_pair(p -> second, p -> first));
		else if(minPQ.top().first < p -> second)
		{
			minPQ.pop();
			minPQ.push(std::make_pair(p -> second, p -> first));
		}
	
	while(!minPQ.empty())
	{
		sampledPairs[minPQ.top().second] = minPQ.top().first;
		minPQ.pop();
	}

	console -> info("Sampled {} color-id pairs. Time-stamp = {}.", sampledPairs.size(),
					time(nullptr) - start_time_);


	auto t_end = time(nullptr);
	console -> info("Sampling abundant color-id pairs took time {} seconds.", t_end - t_start);
	return curBlock;
}



template <typename qf_obj, typename key_obj>
void CdBG_Merger<qf_obj, key_obj>::
	init_disk_buckets()
{
	const uint64_t fileCount1 = cdbg1.get_eq_class_file_count(),
					fileCount2 = cdbg2.get_eq_class_file_count();

	console -> info("Initializing {} x {} disk-buckets.", fileCount1 + 1, fileCount2 + 1);

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
			diskBucket[i][j] = std::ofstream(cdbg.prefix + TEMP_DIR + EQ_ID_PAIRS_FILE + "_" +
											std::to_string(i) + "_" + std::to_string(j)),
			bucketSize[i][j] = 0,
			cumulativeBucketSize[i][j] = 0,
			MPH[i][j] = NULL;
	}

	console -> info("{} x {} disk-buckets  and associated data structures initialized.",
					fileCount1 + 1, fileCount2 + 1);
}



template <typename qf_obj, typename key_obj>
uint64_t CdBG_Merger<qf_obj, key_obj>::
	fill_disk_buckets(uint64_t startingBlock)
{
	auto t_start = time(nullptr);

	console -> info("Writing the non-sampled color-id pairs to disk-files of form ({}). Time-stamp = {}.",
					TEMP_DIR + EQ_ID_PAIRS_FILE + std::string("_X_Y"), time(nullptr) - start_time_);



	uint64_t writtenPairsCount = 0;

	console -> info("Iterating over the CQFs for the non-sampled color-id pairs starting from block {}", startingBlock);

//	std::cerr << "\n\n\nMINIMIZER MAP SIZE= " << minimizerMap.size() << "\n\n\n";
	for (auto & m : minimizerMap) m.clear();

	auto kmerMask = (1ULL << kbits) - 1;

	uint64_t curBlock{startingBlock}, kmerCount{0};
	cdbg1.replaceCQFInMemory(curBlock);
	cdbg2.replaceCQFInMemory(curBlock);
	uint64_t maxMinimizer{0}, minMinimizer{0};
	while(curBlock < cdbg1.get_numBlocks() or curBlock < cdbg2.get_numBlocks()) {
		std::cerr << "\nCurrent Block=" << curBlock << "\n";
		uint64_t maxMinimizer1{0}, minMinimizer1{invalid},
				maxMinimizer2{0}, minMinimizer2{invalid};

		cdbg1.replaceCQFInMemory(curBlock);
		if (curBlock < cdbg1.get_numBlocks()) {
			const CQF<key_obj> *cqf1 = cdbg1.get_current_cqf();
			if (cqf1 == nullptr) {
				std::cerr << "CQF1 is null for block " << curBlock << " while it shouldn't\n";
				std::exit(3);
			}
			typename CQF<key_obj>::Iterator it1 = cqf1->begin();
			uint64_t cntr{0};
			uint64_t count{cqf1->dist_elts()};
			auto minmax = cdbg1.getMinMaxMinimizer(curBlock);
			minMinimizer1 = minmax.first;
			maxMinimizer1 = minmax.second;
			while (!it1.done()) {
				auto keyval = it1.get_cur_hash();
				auto key = hash_64i(keyval.key, kmerMask);
				auto minimizerPair = cdbg1.findMinimizer(key, kbits);
				std::vector<uint64_t> pairs{minimizerPair.first, minimizerPair.second};
				for (auto minimizer : pairs) {
					if (minimizer != invalid and minimizer >= minMinimizer1 and minimizer <= maxMinimizer1) {
						auto &map = minimizerMap[minimizer];
						if (map.find(keyval.key) == map.end()) {
							map[keyval.key] = std::make_pair(keyval.count, 0);
						} else if (map[keyval.key].first and map[keyval.key].first != keyval.count) {
							console->error("found a key with different color class in the second cqf.");
							std::exit(3);
						} else if (not map[keyval.key].first) {
							map[keyval.key].first = keyval.count;
						}
					}
				}
				++it1;
				cntr++;
				if (cntr % PROGRESS_STEP == 0) {
					std::cerr << "\rFirst mantis: " << cntr << " out of " << count;
				}
			}
		}
		cdbg2.replaceCQFInMemory(curBlock);
		if (curBlock < cdbg2.get_numBlocks()) {
			const CQF<key_obj> *cqf2 = cdbg2.get_current_cqf();
			typename CQF<key_obj>::Iterator it2 = cqf2->begin();
			uint64_t cntr{0};
			uint64_t count{cqf2->dist_elts()};
			auto minmax = cdbg2.getMinMaxMinimizer(curBlock);
			minMinimizer2 = minmax.first;
			maxMinimizer2 = minmax.second;
			while (!it2.done()) {
				auto keyval = it2.get_cur_hash();
				auto key = hash_64i(keyval.key, kmerMask);
				auto minimizerPair = cdbg2.findMinimizer(key, kbits);
				std::vector<uint64_t> pairs{minimizerPair.first, minimizerPair.second};
				for (auto minimizer : pairs) {
					if (minimizer != invalid and minimizer >= minMinimizer2 and minimizer <= maxMinimizer2) {
						auto &map = minimizerMap[minimizer];
						if (map.find(keyval.key) == map.end()) {
							map[keyval.key] = std::make_pair(0, keyval.count);
						} else if (map[keyval.key].second and map[keyval.key].second != keyval.count) {
							console->error("found a key with different color class in the second cqf.");
							std::exit(3);
						} else if (not map[keyval.key].second) {
							map[keyval.key].second = keyval.count;
						}
					}
				}
				++it2;
				cntr++;
				if (cntr % PROGRESS_STEP == 0) {
					std::cerr << "\rSecond mantis: " << cntr << " out of " << count;
				}
			}
		}
		minMinimizer = curBlock >= cdbg1.get_numBlocks() ? minMinimizer2 :
				curBlock >= cdbg2.get_numBlocks() ? minMinimizer1 : std::min(minMinimizer1, minMinimizer2);
		maxMinimizer = curBlock >= cdbg1.get_numBlocks() ? maxMinimizer2 :
					   curBlock >= cdbg2.get_numBlocks() ? maxMinimizer1 : std::min(maxMinimizer1, maxMinimizer2);
		std::cerr << "\nMin minimizer=" << minMinimizer << " Max minimizer=" << maxMinimizer << "\n";
		curBlock++;
		uint64_t b = minMinimizer;
		while (b <= maxMinimizer) {
			if (not minimizerMap[b].empty()) {
				for (auto &kv : minimizerMap[b]) {
					colorCount1 = std::max(colorCount1, static_cast<uint64_t >(kv.second.first)),
					colorCount2 = std::max(colorCount2, static_cast<uint64_t >(kv.second.second));

					if(sampledPairs.find(kv.second) == sampledPairs.end())
					{
						add_color_id_pair(kv.second.first, kv.second.second, diskBucket);
						writtenPairsCount++;
					}
					kmerCount++;
					if(kmerCount % PROGRESS_STEP == 0)
						std::cerr << "\rObserved count of distinct k-mers: " << kmerCount * 10 / PROGRESS_STEP
						<< "M, written color-id pairs to disk: " << writtenPairsCount
						<< ". Time-stamp = " << time(nullptr) - start_time_ << ".";
				}
			}
			minimizerMap[b].clear();
			b++;
		}
		std::cerr << "\r";
	}
	console -> info("Distinct kmers found {}, color-id pairs written to disk {}. Time-stamp = {}.",
					kmerCount, writtenPairsCount, time(nullptr) - start_time_);


	// Flush and close the disk-bucket files.
	for(int i = 0; i <= cdbg1.get_eq_class_file_count(); ++i)
		for(int j = 0; j <= cdbg2.get_eq_class_file_count(); ++j)
		{
			diskBucket[i][j].flush();
			diskBucket[i][j].close();
		}

	
	auto t_end = time(nullptr);
	console -> info("Filling up the disk-buckets with color-id pairs took time {} seconds.", t_end - t_start);


	return kmerCount;
}



template <typename qf_obj, typename key_obj>
void CdBG_Merger<qf_obj, key_obj>::
	add_color_id_pair(const uint64_t colorID1, const uint64_t colorID2,
						std::vector<std::vector<std::ofstream>> &diskBucket)
{
	// TODO: Add faster file-write mechanism.

	const uint64_t row = (colorID1 ? (colorID1 - 1) / mantis::NUM_BV_BUFFER + 1 : 0),
					col = (colorID2 ? (colorID2 - 1) / mantis::NUM_BV_BUFFER + 1 : 0);

	diskBucket[row][col] << colorID1 << " " << colorID2 << "\n";
}



template <typename qf_obj, typename key_obj>
uint64_t CdBG_Merger<qf_obj, key_obj>::
	filter_disk_buckets()
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

	uint maxMemoryForSort = std::max(get_max_sort_memory(), (uint64_t)1);
	
	for(int i = 0; i <= fileCount1; ++i)
		for(int j = 0; j <= fileCount2; ++j)
		{
			std::string diskBucket = cdbg.prefix + TEMP_DIR + EQ_ID_PAIRS_FILE +
									"_" + std::to_string(i) + "_" + std::to_string(j);

			std::string sysCommand = "sort -u";
			sysCommand += " --parallel=" + std::to_string(threadCount);
			sysCommand += " -S " + std::to_string(maxMemoryForSort) + "G";
			sysCommand += " -o " + diskBucket + " " + diskBucket;

			console -> info("System command used:\n{}", sysCommand);

			system(sysCommand.c_str());

			
			std::string lineCountFile = cdbg.prefix + TEMP_DIR + ID_PAIR_COUNT_FILE;
			sysCommand = "wc -l " + diskBucket + " | egrep -o \"[0-9]+ \" > " + lineCountFile;
			system(sysCommand.c_str());

			std::ifstream lineCount(lineCountFile);

			lineCount >> bucketSize[i][j];
			lineCount.close();

			colorClassCount += bucketSize[i][j];

			console -> info("Filtered {} unique color-id pairs at disk-bucket {}. Time-stamp = {}.",
							bucketSize[i][j], diskBucket, time(nullptr) - start_time_);
		}


	console -> info("Count of unique color-id pairs = {}. Time-stamp = {}.",
					colorClassCount, time(nullptr) - start_time_);


	auto t_end = time(nullptr);
	console -> info("Filtering the unique color-id pairs took time {} seconds.", t_end - t_start);

	return colorClassCount;
}



template <typename qf_obj, typename key_obj>
uint64_t CdBG_Merger<qf_obj, key_obj>::
	get_max_sort_memory()
{
	uint64_t bvBuffMemory = mantis::NUM_BV_BUFFER * cdbg.num_samples / 8;
	uint64_t maxRRR1size = 0, maxRRR2size = 0;

	// File-size calculation reference:
	// https://stackoverflow.com/questions/5840148/how-can-i-get-a-files-size-in-c

	// Determine the maximum bitvectorRRR file-size for cdbg1.
	for(uint64_t i = 0; i < cdbg1.get_eq_class_file_count(); ++i)
	{
		struct stat64 stat_buf;
		if(stat64(cdbg1.get_eq_class_files()[i].c_str(), &stat_buf) == 0)
			maxRRR1size = std::max(maxRRR1size, (uint64_t)stat_buf.st_size);
		else
		{
			console -> error("File size of the bitvectorRRR file {} for CdBG1 cannot be determined.",
							cdbg1.get_eq_class_files()[i]);
			exit(1);
		}
	}

	// Determine the maximum bitvectorRRR file-size for cdbg2.
	for(uint64_t i = 0; i < cdbg2.get_eq_class_file_count(); ++i)
	{
		struct stat64 stat_buf;
		if(stat64(cdbg2.get_eq_class_files()[i].c_str(), &stat_buf) == 0)
			maxRRR2size = std::max(maxRRR2size, (uint64_t)stat_buf.st_size);
		else
		{
			console -> error("File size of the bitvectorRRR file {} for CdBG2 cannot be determined.",
							cdbg2.get_eq_class_files()[i]);
			exit(1);
		}
	}


	return (bvBuffMemory + maxRRR1size + maxRRR2size) / (1024 * 1024 * 1024);
}



template <typename qf_obj, typename key_obj>
void CdBG_Merger<qf_obj, key_obj>:: assign_abundant_color_ids()
{
	auto t_start = time(nullptr);

	console -> info("Assigning color-ids to the {} sampled (on abundance) color-id pairs. Time-stamp = {}.",
					sampledPairs.size(), time(nullptr) - start_time_);


	std::vector<std::pair<uint64_t, std::pair<uint64_t, uint64_t>>> idPairs;
	idPairs.reserve(sampledPairs.size());

	for(auto it = sampledPairs.begin(); it != sampledPairs.end(); ++it)
		idPairs.emplace_back(it -> second, it -> first);

	sort(idPairs.begin(), idPairs.end(), std::greater<std::pair<uint64_t, std::pair<uint64_t, uint64_t>>>());
	
	uint64_t colorId = 0;
	for(auto it = idPairs.begin(); it != idPairs.end(); ++it)
		sampledPairs[it -> second] = ++colorId;


	auto t_end = time(nullptr);
	console -> info("Assigning color-ids to the sampled pairs. Took time {} seconds.", t_end - t_start);
}



template <typename qf_obj, typename key_obj>
void CdBG_Merger<qf_obj, key_obj>::
	build_MPH_tables()
{	
	auto t_start = time(nullptr);

	console -> info("Building an MPH (Minimal Perfect Hash) table per disk-bucket. Time-stamp = {}.",
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


				console -> info("Loading the unique color-id pairs from bucket ({}, {}) into memory. Time-stamp = {}.",
								i, j, time(nullptr) - start_time_);

				std::ifstream input(cdbg.prefix + TEMP_DIR + EQ_ID_PAIRS_FILE +
									"_" + std::to_string(i) + "_" + std::to_string(j));
				uint64_t id1, id2;

				// TODO: Add faster file-read mechanism.
				while(input >> id1 >> id2)
					colorIdPair.emplace_back(id1, id2);

				input.close();

				console -> info("Loaded {} color-id pairs into memory. Time-stamp = {}.", colorIdPair.size(),
								time(nullptr) - start_time_);

				
				console -> info("Constructing a BooPHF with {} elements using {} threads. Time-stamp = {}.",
								colorIdPair.size(), threadCount, time(nullptr) - start_time_);

				// Lowest bit/elem is achieved with gamma=1, higher values lead to larger mphf but faster construction/query.
				double gammaFactor = 2.0;	// gamma = 2 is a good tradeoff (leads to approx 3.7 bits/key).

				// Build the mphf.
				MPH[i][j] = new boophf_t(colorIdPair.size(), colorIdPair, threadCount, gammaFactor);
				mphMemory += MPH[i][j] -> totalBitSize();

				
				console -> info("For bucket ({}, {}) BooPHF constructed perfect hash for {} keys; total memory = {} MB, bits/elem : {}. Time-stamp = {}.",
								i, j, colorIdPair.size(), (MPH[i][j] -> totalBitSize() / 8) / (1024 * 1024),
								(double)(MPH[i][j] -> totalBitSize()) / colorIdPair.size(),
                                time(nullptr) - start_time_);

				
				colorIdPair.clear();
				colorIdPair.shrink_to_fit();
			}

	console -> info("Total memory consumed by all the MPH tables = {} MB.", (mphMemory / 8) / (1024 * 1024));
	

	auto t_end = time(nullptr);
	console -> info("Building the MPH tables took time {} seconds.", t_end - t_start);
}



template <typename qf_obj, typename key_obj>
void CdBG_Merger<qf_obj, key_obj> ::
	concat(const BitVectorRRR &bv1, const uint64_t colCount1, const uint64_t colorID1,
			const BitVectorRRR &bv2, const uint64_t colCount2, const uint64_t colorID2,
			BitVector &resultVec)
{
	const uint64_t wordLen = 64;


	if(colorID1)	// Color ID = 0 implies an absent color-id.
	{
		uint64_t offset = ((colorID1 - 1) % mantis::NUM_BV_BUFFER) * colCount1;

		for(uint32_t blockStart = 0; blockStart < (colCount1 / wordLen) * wordLen; blockStart += wordLen)
				resultVec.set_int(blockStart, bv1.get_int(offset + blockStart, wordLen), wordLen);

		if(colCount1 % wordLen)
			resultVec.set_int((colCount1 / wordLen) * wordLen,
								bv1.get_int(offset + ((colCount1 / wordLen) * wordLen), colCount1 % wordLen),
								colCount1 % wordLen);
	}


	if(colorID2)	// Color ID = 0 implies an absent color-id.
	{
		uint64_t offset = ((colorID2 - 1) % mantis::NUM_BV_BUFFER) * colCount2;

		for(uint32_t blockStart = 0; blockStart < (colCount2 / wordLen) * wordLen; blockStart += wordLen)
				resultVec.set_int(colCount1 + blockStart, bv2.get_int(offset + blockStart, wordLen), wordLen);

		if(colCount2 % wordLen)
			resultVec.set_int(colCount1 + (colCount2 / wordLen) * wordLen,
								bv2.get_int(offset + ((colCount2 / wordLen) * wordLen), colCount2 % wordLen),
								colCount2 % wordLen);
	}
}



template <typename qf_obj, typename key_obj>
void CdBG_Merger<qf_obj, key_obj>::
	build_color_class_table()
{
	auto t_start = time(nullptr);

	console -> info("At color-class building (bitvectors concatenation) phase. Time-stamp = {}.",
					time(nullptr) - start_time_);

	
	uint64_t writtenPairsCount = 0;

	console -> info("Building color-class for the {} sampled (on abundance) color-id pairs. Time-stamp = {}.",
					sampledPairs.size(), time(nullptr) - start_time_);

	build_abundant_color_classes();
	
	writtenPairsCount = sampledPairs.size();
	if(writtenPairsCount % mantis::NUM_BV_BUFFER == 0)
		bv_buffer_serialize(writtenPairsCount);

	console -> info("Color-classes built for the {} sampled (on abundance) color-id pairs. Time-stamp = {}.",
					sampledPairs.size(), time(nullptr) - start_time_);


	const uint64_t fileCount1 = cdbg1.get_eq_class_file_count(), fileCount2 = cdbg2.get_eq_class_file_count();
	BitVectorRRR bitVec1, bitVec2;


	for(uint64_t i = 0; i <= fileCount1; ++i)
	{
		if(i > 0)
		{
			sdsl::load_from_file(bitVec1, cdbg1.get_eq_class_files()[i - 1]);
			console -> info("Mantis 1: loaded one bitvectorRRR from file {}. Time-stamp = {}.",
							cdbg1.get_eq_class_files()[i - 1], time(nullptr) - start_time_);
		}

		for(uint64_t j = 0; j <= fileCount2; ++j)
		{
			if(j > 0)
			{
				sdsl::load_from_file(bitVec2, cdbg2.get_eq_class_files()[j - 1]);
				console -> info("Mantis 2: loaded one bitvectorRRR from file {}. Time-stamp = {}.",
								cdbg2.get_eq_class_files()[j - 1], time(nullptr) - start_time_);
			}

			console -> info("At bucket ({}, {}), size = {}. Time-stamp = {}.", i, j, bucketSize[i][j],
							time(nullptr) - start_time_);


			cumulativeBucketSize[i][j] = writtenPairsCount;


			if(bucketSize[i][j])
			{
				uint64_t queueCount = ((cumulativeBucketSize[i][j] + bucketSize[i][j] - 1) / mantis::NUM_BV_BUFFER)
										- (cumulativeBucketSize[i][j] / mantis::NUM_BV_BUFFER) + 1;

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


				std::ifstream input(cdbg.prefix + TEMP_DIR + EQ_ID_PAIRS_FILE +
									"_" + std::to_string(i) + "_" + std::to_string(j));
				std::pair<uint64_t, uint64_t> idPair;

				// TODO: Add faster file-read mechanism.
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
						BitVector colorClass(cdbg.num_samples);

						concat(bitVec1, cdbg1.get_num_samples(), it -> first,
								bitVec2, cdbg2.get_num_samples(), it -> second, colorClass);

						cdbg.add_bitvector(colorClass, colorID - 1),
						writtenPairsCount++;

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


	// Serialize the bitvector buffer last time if needed.
	if (writtenPairsCount % mantis::NUM_BV_BUFFER > 0)
	{
		console -> info("Serializing bitvector buffer with {} color-classes.", writtenPairsCount);
		bv_buffer_serialize(writtenPairsCount);
	}


	auto t_end = time(nullptr);
	console -> info("Color-class building phase took time {} seconds.", t_end - t_start);
}



template <typename qf_obj, typename key_obj>
void CdBG_Merger<qf_obj, key_obj>::
	build_abundant_color_classes()
{
	const uint64_t fileCount1 = cdbg1.get_eq_class_file_count(), fileCount2 = cdbg2.get_eq_class_file_count();
	BitVectorRRR bitVec1, bitVec2;

	std::vector<std::pair<uint64_t, uint64_t>> idPairs;
	
	idPairs.reserve(sampledPairs.size());
	for(auto it = sampledPairs.begin(); it != sampledPairs.end(); ++it) {
		idPairs.emplace_back(it->first);
	}

	sort(idPairs.begin(), idPairs.end(),
		[](const std::pair<uint64_t, uint64_t> &lhs, const std::pair<uint64_t, uint64_t> &rhs)
		{
			if(lhs.first != rhs.first)
				return lhs.first < rhs.first;
			return lhs.second < rhs.second;
		});

	/*for (auto fs : idPairs) {
		std::cerr << fs.first << " " << fs.second << "\n";
	}*/
	uint64_t currBucket1{std::numeric_limits<uint64_t>::max()}, currBucket2{std::numeric_limits<uint64_t>::max()};

	for(auto fs : idPairs)
	{
		uint64_t reqBucket1 = fs.first / mantis::NUM_BV_BUFFER,
					reqBucket2 = fs.second / mantis::NUM_BV_BUFFER;

		if(/*reqBucket1 && */reqBucket1 != currBucket1)
		{
			sdsl::load_from_file(bitVec1, cdbg1.get_eq_class_files()[reqBucket1]);
			currBucket1 = reqBucket1;

			console -> info("Mantis 1: loaded one bitvectorRRR from file {}. Time-stamp = {}.",
							cdbg1.get_eq_class_files()[currBucket1], time(nullptr) - start_time_);
		}

		if(/*reqBucket2 && */reqBucket2 != currBucket2)
		{
			sdsl::load_from_file(bitVec2, cdbg2.get_eq_class_files()[reqBucket2]);
			currBucket2 = reqBucket2;

			console -> info("Mantis 2: loaded one bitvectorRRR from file {}. Time-stamp = {}.",
							cdbg2.get_eq_class_files()[currBucket2], time(nullptr) - start_time_);
		}


		uint64_t colorID = sampledPairs[fs];
		BitVector colorClass(cdbg.num_samples);
		concat(bitVec1, cdbg1.get_num_samples(), fs.first,
				bitVec2, cdbg2.get_num_samples(), fs.second, colorClass);
//		std::cerr << colorID << "\n";
		cdbg.add_bitvector(colorClass, colorID-1);
		/*
			No bitvector-buffer serialization required, as the number of sampled pairs'
			count (SAMPLE_PAIR_COUNT) is defined as less than or equal to the bitvector-buffer size.
			(Check its declaration)
		*/
	}
}



template <class qf_obj, class key_obj>
void CdBG_Merger<qf_obj, key_obj>:: bv_buffer_serialize(uint64_t colorClsCount)
{
	if (colorClsCount % mantis::NUM_BV_BUFFER > 0)
		cdbg.bv_buffer.resize((colorClsCount % mantis::NUM_BV_BUFFER) * cdbg.num_samples);
	
	BitVectorRRR final_com_bv(cdbg.bv_buffer);
	std::string bv_file(cdbg.prefix + std::to_string(cdbg.num_serializations) + "_" + mantis::EQCLASS_FILE);
	sdsl::store_to_file(final_com_bv, bv_file);

	cdbg.num_serializations++;
}



template <typename qf_obj, typename key_obj>
void CdBG_Merger<qf_obj, key_obj> ::
	initialize_CQF(uint32_t keybits, qf_hashmode hashMode, uint32_t seed, uint64_t kmerCount)
{
	// Get floor(log2(kmerCount))
	uint32_t qbits;
	for(qbits = 0; (kmerCount >> qbits) != (uint64_t)1; qbits++);

	// Get ceil(log2(kmerCount))
	if(kmerCount & (kmerCount - 1))	// if kmerCount is not a power of 2
		qbits++;
	
	qbits += 2;	// to avoid the initial rapid resizes at minuscule load factors

	
	if(cdbg.dbg_alloc_flag == MANTIS_DBG_IN_MEMORY)
	{
		CQF<key_obj> cqf(qbits, keybits, hashMode, seed);
		cdbg.dbg = cqf;
	}
	else if(cdbg.dbg_alloc_flag == MANTIS_DBG_ON_DISK)
	{
		CQF<key_obj> cqf(qbits, keybits, hashMode, seed, cdbg.prefix + mantis::CQF_FILE);
		cdbg.dbg = cqf;
	}
	else
	{
		ERROR("Wrong Mantis alloc mode.");
		exit(EXIT_FAILURE);
	}

	cdbg.dbg.set_auto_resize();
}


template <typename qf_obj, typename key_obj>
void CdBG_Merger<qf_obj, key_obj>::build_CQF()
{
	auto t_start = time(nullptr);

	console -> info("At CQFs merging phase. Time-stamp = {}.\n", time(nullptr) - start_time_);


	uint64_t kmerCount{0}, foundAbundantId{0};
	cdbg1.replaceCQFInMemory(0);
	cdbg2.replaceCQFInMemory(0);


	auto kmerMask = (1ULL << kbits) - 1;

	uint64_t curBlock{0}, outputCQFBlockId{0}, blockKmerCnt{0};
	for (auto &m : cdbg.minimizerCntr) m = 0;

	auto hash_mode = cdbg1.get_current_cqf()->hash_mode();
	auto seed = cdbg1.get_current_cqf()->seed();
	cdbg.initializeNewCQFBlock(outputCQFBlockId, kbits, hash_mode, seed);
	uint64_t maxMinimizer{0}, minMinimizer{invalid};
	while(curBlock < cdbg1.get_numBlocks() or curBlock < cdbg2.get_numBlocks()) {
		std::cerr << "Current block=" << curBlock << "\n";
		uint64_t maxMinimizer1{0}, minMinimizer1{invalid},
				maxMinimizer2{0}, minMinimizer2{invalid};

		cdbg1.replaceCQFInMemory(curBlock);
		if (curBlock < cdbg1.get_numBlocks()) {
			const CQF<key_obj> *cqf1 = cdbg1.get_current_cqf();
			typename CQF<key_obj>::Iterator it1 = cqf1->begin();
            uint64_t cntr{0}, count{cqf1->dist_elts()};
			auto minmax = cdbg1.getMinMaxMinimizer(curBlock);
			minMinimizer1 = minmax.first;
			maxMinimizer1 = minmax.second;
			while (!it1.done()) {
				auto keyval = it1.get_cur_hash();
				auto key = hash_64i(keyval.key, kmerMask);
				auto minimizerPair = cdbg1.findMinimizer(key, kbits);
				std::vector<uint64_t> pairs{minimizerPair.first, minimizerPair.second};
				for (auto minimizer : pairs) {
					if (minimizer != invalid and minimizer >= minMinimizer1 and minimizer <= maxMinimizer1) {
						auto &map = minimizerMap[minimizer];
						if (map.find(keyval.key) == map.end()) {
							cdbg.minimizerCntr[minimizer]++;
							map[keyval.key] = std::make_pair(keyval.count, 0);
						} else if (map[keyval.key].first and map[keyval.key].first != keyval.count) {
							console->error("found a key with different color class in the second cqf.");
							std::exit(3);
						} else if (not map[keyval.key].first) {
							map[keyval.key].first = keyval.count;
						}
					}
				}
				++it1;
                cntr++;
                if (cntr % 10000000 == 0) {
                    std::cerr << "\r" << cntr << " out of " << count;
                }
			}
		}

		cdbg2.replaceCQFInMemory(curBlock);
		if (curBlock < cdbg2.get_numBlocks()) {
			const CQF<key_obj> *cqf2 = cdbg2.get_current_cqf();
			typename CQF<key_obj>::Iterator it2 = cqf2->begin();
            uint64_t cntr{0}, count{cqf2->dist_elts()};
			auto minmax = cdbg2.getMinMaxMinimizer(curBlock);
			minMinimizer2 = minmax.first;
			maxMinimizer2 = minmax.second;
			while (!it2.done()) {
				auto keyval = it2.get_cur_hash();
				auto key = hash_64i(keyval.key, kmerMask);
				auto minimizerPair = cdbg2.findMinimizer(key, kbits);
				std::vector<uint64_t> pairs{minimizerPair.first, minimizerPair.second};
				for (auto minimizer : pairs) {
					if (minimizer != invalid and minimizer >= minMinimizer2 and minimizer <= maxMinimizer2) {
						auto &map = minimizerMap[minimizer];
						if (map.find(keyval.key) == map.end()) {
							cdbg.minimizerCntr[minimizer]++;
							map[keyval.key] = std::make_pair(0, keyval.count);
						} else if (map[keyval.key].second and map[keyval.key].second != keyval.count) {
							console->error("found a key with different color class in the second cqf.");
							std::exit(3);
						} else if (not map[keyval.key].second) {
							map[keyval.key].second = keyval.count;
						}
					}
				}
				++it2;
                cntr++;
                if (cntr % 10000000 == 0) {
                    std::cerr << "\r" << cntr << " out of " << count;
                }
            }
		}
		minMinimizer = curBlock >= cdbg1.get_numBlocks() ? minMinimizer2 :
					   curBlock >= cdbg2.get_numBlocks() ? minMinimizer1 : std::min(minMinimizer1, minMinimizer2);
		maxMinimizer = curBlock >= cdbg1.get_numBlocks() ? maxMinimizer2 :
					   curBlock >= cdbg2.get_numBlocks() ? maxMinimizer1 : std::min(maxMinimizer1, maxMinimizer2);
		std::cerr << "\nMin minimizer=" << minMinimizer << " Max minimizer=" << maxMinimizer << "\n";
		curBlock++;

        uint64_t m = minMinimizer;
        while (m <= maxMinimizer) {
            if (not minimizerMap[m].empty()) {
				if (m > 0 and blockKmerCnt + minimizerMap[m].size() > block_kmer_threshold) {
					blockKmerCnt = 0;
					cdbg.serializeCurrentCQF();
					outputCQFBlockId++;
					cdbg.initializeNewCQFBlock(outputCQFBlockId, kbits, hash_mode, seed);
				}
				blockKmerCnt += minimizerMap[m].size();
				for (auto &kv : minimizerMap[m]) {
                    uint64_t colorId = get_color_id(kv.second);
					auto keyObj = KeyObject(kv.first, 0, colorId);
					cdbg.add_kmer2CurDbg(keyObj, QF_NO_LOCK | QF_KEY_IS_HASH);
                    kmerCount++;

                    if(colorId <= sampledPairs.size())
                        foundAbundantId++;

                    if(kmerCount % PROGRESS_STEP == 0)
						std::cerr << "\rKmers merged: " << kmerCount * 10 / PROGRESS_STEP
								  << "M, Time-stamp = " << time(nullptr) - start_time_ << "      ";
                }
            }
            cdbg.minimizerBorder[m] = outputCQFBlockId;
            minimizerMap[m].clear();
            m++;
        }
        std::cerr << "\r";
	}
	cdbg.serializeCurrentCQF();

	console -> info("Total kmers merged: {}. Time-stamp: {}.", kmerCount, time(nullptr) - start_time_);

	console -> info("Out of {} kmers, {} have color-ids that are sampled earlier.", kmerCount, foundAbundantId);


	auto t_end = time(nullptr);
	console -> info("Merging the CQFs took time {} seconds.", t_end - t_start);
}



template <typename qf_obj, typename key_obj>
uint64_t CdBG_Merger<qf_obj, key_obj>:: get_color_id(const std::pair<uint64_t, uint64_t> &idPair)
{
	auto it = sampledPairs.find(idPair);
	if(it != sampledPairs.end())
		return it -> second;

	const uint64_t row = (idPair.first ? (idPair.first - 1) / mantis::NUM_BV_BUFFER + 1 : 0),
					col = (idPair.second ? (idPair.second - 1) / mantis::NUM_BV_BUFFER + 1 : 0);


	return cumulativeBucketSize[row][col] + MPH[row][col] -> lookup(idPair) + 1;
}

template <typename qf_obj, typename key_obj>
inline void CdBG_Merger<qf_obj, key_obj> ::
	add_kmer(uint64_t kmer, uint64_t colorID, uint64_t &step, typename CQF<key_obj>::Iterator &walkBehindIterator)
{
	if(cdbg.dbg.insert(KeyObject(kmer, 0, colorID), QF_NO_LOCK | QF_KEY_IS_HASH) == QF_NO_SPACE)
	{
		// Auto_resize failed.
		console -> error("The CQF is full and auto resize failed. Please re-run build with a bigger size.");
		exit(1);
	}

	
	step++;
	if(step == ITERATOR_WINDOW_SIZE)
		walkBehindIterator = cdbg.dbg.begin(true);
	else if(step > ITERATOR_WINDOW_SIZE)
		++walkBehindIterator;
}



template <typename qf_obj, typename key_obj>
void CdBG_Merger<qf_obj, key_obj>:: serialize_cqf_and_sampleid_list()
{
	cdbg.serializeCurrentCQF();
	// Serialize the CQF.
/*
	if(cdbg.dbg_alloc_flag == MANTIS_DBG_IN_MEMORY)
		cdbg.dbg.serialize(cdbg.prefix + mantis::CQF_FILE);
	else
		cdbg.dbg.close();
*/

	console -> info("Serialized CQF.");


	// Serialize the sample-id map.
	std::ofstream outputFile(cdbg.prefix + mantis::SAMPLEID_FILE);

	for(auto idSample : cdbg.sampleid_map)
		outputFile << idSample.first << " " << idSample.second << "\n";

	outputFile.close();

	console -> info("Serialized sample-id mapping.");

	std::ofstream minfile(cdbg.prefix + mantis::MINIMIZER_FREQ, std::ios::binary);
	minfile.write(reinterpret_cast<char *>(cdbg.minimizerCntr.data()), cdbg.minimizerCntr.size()*sizeof(typename decltype(cdbg.minimizerCntr)::value_type));
	minfile.close();

	minfile.open(cdbg.prefix + mantis::MINIMIZER_BOUNDARY, std::ios::binary);
	minfile.write(reinterpret_cast<char *>(cdbg.minimizerBorder.data()),
				  cdbg.minimizerBorder.size()*sizeof(typename decltype(cdbg.minimizerBorder)::value_type));
	minfile.close();
}

template <typename qf_obj, typename key_obj>
void CdBG_Merger<qf_obj, key_obj>::
calc_mst_stats(ColoredDbg<qf_obj, key_obj> &cdbg1, ColoredDbg<qf_obj, key_obj> &cdbg2, std::string& dir1, std::string& dir2)
{
	auto t_start = time(nullptr);

	console -> info("At color-class building (bitvectors concatenation) phase. Time-stamp = {}.",
					time(nullptr) - start_time_);


	const uint64_t bucketCnt1 = cdbg1.get_eq_class_file_count(), bucketCnt2 = cdbg2.get_eq_class_file_count();

	uint64_t writtenPairsCount = 0;
	BitVectorRRR bitVec1, bitVec2;

	std::cerr << dir1 + mantis::PARENTBV_FILE << "\n";
	std::cerr << dir2 + mantis::PARENTBV_FILE << "\n";
	sdsl::int_vector<> parentbv1;
	sdsl::load_from_file(parentbv1, dir1 + "/" + mantis::PARENTBV_FILE);
	sdsl::int_vector<> parentbv2;
	sdsl::load_from_file(parentbv2, dir2 + "/" + mantis::PARENTBV_FILE);

	std::cerr << "loaded the files\n";
	spp::sparse_hash_map<std::pair<uint64_t, uint64_t >, bool, Custom_Pair_Hasher> isVirtual;

	std::cerr << "parents sizes: " << parentbv1.size() << " " << parentbv2.size() << "\n";
	uint64_t zero1 = parentbv1.size() - 1;
	uint64_t zero2 = parentbv2.size() - 1;

	std::cerr << "# of buckets: " << bucketCnt1 << " " << bucketCnt2 << "\n";
	for (uint64_t i = 0; i <= bucketCnt1; ++i)
	{
		std::cerr << "\n" << i;
		for(uint64_t j = 0; j <= bucketCnt2; ++j) {
			std::cerr << "\n" << j << "\n";
			std::ifstream input(cdbg.prefix + TEMP_DIR + EQ_ID_PAIRS_FILE +
								"_" + std::to_string(i) + "_" + std::to_string(j));
			std::pair<uint64_t, uint64_t> idPair;
			uint64_t cntr{0};
			while (input >> idPair.first >> idPair.second) {
//				std::cerr << "reading the file: " << idPair.first << " " << idPair.second << "\n";
				if (!idPair.first) idPair.first = zero1;
				else idPair.first--;
				if (!idPair.second) idPair.second = zero2;
				else idPair.second--;
				// doesn't matter if the idPair is in the map or not, we want to set the value to false
				isVirtual[idPair] = false;
				uint64_t pid1{idPair.first}, cid1{idPair.first}, pid2{idPair.second}, cid2{idPair.second};
				std::cerr << "\r" << cntr++;
				do {
//					std::cerr << cid1 << " " << cid2 << "\n";
					cid1 = pid1;
					cid2 = pid2;
					pid1 = parentbv1[cid1];
					pid2 = parentbv2[cid2];
					auto ppair = std::make_pair(pid1, pid2);
					if (isVirtual.find(ppair) == isVirtual.end()) {
						isVirtual[ppair] = true;
					}
				} while (pid1 != cid1 and pid2 != cid2);
			}
		}
	}

	uint64_t actualCnt{0}, virtualCnt{0};
	for (auto& kv : isVirtual) {
		if (kv.second) virtualCnt++;
		else actualCnt++;
	}

	std::cerr << "actual cnt: " << actualCnt << " virtual cnt: " << virtualCnt << "\n";
	std::exit(1);
}

template <typename qf_obj, typename key_obj>
uint64_t CdBG_Merger<qf_obj, key_obj>::
store_abundant_color_pairs(std::ofstream& output)
{
	for (auto &idpair : sampledPairs) {
		uint64_t colorID = idpair.second-1;
		auto fs = idpair.first;
		output.write(reinterpret_cast<char*>(&colorID), sizeof(colorID));
		output.write(reinterpret_cast<char*>(&(fs.first)), sizeof(fs.first));
		output.write(reinterpret_cast<char*>(&(fs.second)), sizeof(fs.second));
//		std::cerr << colorID << " " << fs.first << " " << fs.second << "\n";
	}
	return sampledPairs.size();
}


template <typename qf_obj, typename key_obj>
void CdBG_Merger<qf_obj, key_obj>::
store_color_pairs(ColoredDbg<qf_obj, key_obj> &cdbg1, ColoredDbg<qf_obj, key_obj> &cdbg2,
				  uint64_t& numColorBuffers)
{
	auto t_start = time(nullptr);

	console -> info("At color-class building (bitvectors concatenation) phase. Time-stamp = {}.",
					time(nullptr) - start_time_);


	const uint64_t fileCount1 = cdbg1.get_eq_class_file_count(), fileCount2 = cdbg2.get_eq_class_file_count();

	console->info("EQClass Count for cdbg1 and 2: {}, {}", fileCount1, fileCount2);
	uint64_t writtenPairsCount = 0;
	//BitVectorRRR bitVec1, bitVec2;

	console->info("Writing color pairs and ids into file {}", cdbg.prefix + "/newID2oldIDs" );
	std::ofstream output(cdbg.prefix + "newID2oldIDs");
	output.write(reinterpret_cast<char*>(&writtenPairsCount), sizeof(writtenPairsCount));
	writtenPairsCount = store_abundant_color_pairs(output);

	for(uint64_t i = 0; i <= fileCount1; ++i) {

		for (uint64_t j = 0; j <= fileCount2; ++j) {

			console->info("At bucket ({}, {}), size = {}. Time-stamp = {}.", i, j, bucketSize[i][j],
						  time(nullptr) - start_time_);

			std::ifstream input(cdbg.prefix + TEMP_DIR + EQ_ID_PAIRS_FILE +
								"_" + std::to_string(i) + "_" + std::to_string(j));
			std::pair<uint64_t, uint64_t> idPair;
			cumulativeBucketSize[i][j] = writtenPairsCount;
			while (input >> idPair.first >> idPair.second) {
				/// @fatemeh write down the pair and the associated colorID here
				uint64_t colorID = cumulativeBucketSize[i][j] + MPH[i][j]->lookup(idPair);// + 1;
				output.write(reinterpret_cast<char*>(&colorID), sizeof(colorID));
				output.write(reinterpret_cast<char*>(&(idPair.first)), sizeof(idPair.first));
				output.write(reinterpret_cast<char*>(&(idPair.second)), sizeof(idPair.second));
//				std::cerr << colorID << " " << idPair.first << " " << idPair.second << "\n";

				writtenPairsCount++;
			}
			input.close();
		}
	}
	output.seekp(0, std::ios::beg);
	output.write(reinterpret_cast<char*>(&writtenPairsCount), sizeof(writtenPairsCount));
	output.close();
	numColorBuffers = writtenPairsCount/mantis::NUM_BV_BUFFER + 1;

	auto t_end = time(nullptr);
	console -> info("Writing {} color pairs took time {} seconds.", writtenPairsCount, t_end - t_start);
}



template <typename qf_obj, typename key_obj>
void CdBG_Merger<qf_obj, key_obj>::merge()
{

    auto t_start = time(nullptr);
    console -> info ("Splitting output minimizers into blocks based on sum of the two input minimizers");
//    divide_minimizers_into_blocks(cdbg1.minimizerCntr, cdbg2.minimizerCntr);
	minimizerMap.resize(cdbg1.minimizerBorder.size());
	cdbg.minimizerCntr.resize(cdbg1.minimizerCntr.size(), 0);
	cdbg.minimizerBorder.resize(cdbg1.minimizerBorder.size(), 0);
	for (auto &m : minimizerMap) m.clear();
    console -> info ("Done dividing. Time-stamp = {}.\n", time(nullptr) - t_start);

    t_start = time(nullptr);
	console -> info ("Merge starting. Time-stamp = {}.\n", time(nullptr) - start_time_);


	// Make the temporary directory if it doesn't exist.
	std::string tempDir = cdbg.prefix + TEMP_DIR;

	if(!mantis::fs::DirExists(tempDir.c_str()))
		mantis::fs::MakeDir(tempDir.c_str());

	// Check to see if the temporary directory exists now.
	if(!mantis::fs::DirExists(tempDir.c_str()))
	{
		console -> error("Temporary directory {} could not be created.", tempDir);
		exit(1);
	}


	auto tillBlock = sample_color_id_pairs(mantis::SAMPLE_SIZE);

	init_disk_buckets();
	uint64_t kmerCount = fill_disk_buckets(tillBlock);
	uint64_t colorClassCount = sampledPairs.size() + filter_disk_buckets();

	assign_abundant_color_ids();
	build_MPH_tables();

	uint64_t num_colorBuffers = 1;
//	cdbg.bv_buffer = BitVector(mantis::NUM_BV_BUFFER * cdbg.num_samples);
//	build_color_class_table();
//	cdbg.bv_buffer = BitVector(0);

	//	calc_mst_stats(cdbg1, cdbg2, opt.dir1, opt.dir2);
	store_color_pairs(cdbg1, cdbg2, num_colorBuffers);

//	initialize_CQF(cdbg1.get_cqf() -> keybits(), cdbg1.get_cqf() -> hash_mode(), cdbg1.get_cqf() -> seed(), kmerCount);
	build_CQF();


	// Remove the temporary directory.
	std::string sysCommand = "rm -rf " + tempDir;
	console -> info("Removing the temporary directory. System command used:\n{}", sysCommand);

	system(sysCommand.c_str());


	auto t_end = time(nullptr);

	console -> info("Merge completed.");

	console -> info("Input colored dBG 1: over {} samples and has {} k-mers and {} color-classes.",
					cdbg1.get_num_samples(), cdbg1.dbg.dist_elts(), colorCount1);
	console -> info("Input colored dBG 2: over {} samples and has {} k-mers and {} color-classes.",
					cdbg2.get_num_samples(), cdbg2.dbg.dist_elts(), colorCount2);
	console -> info("Merged colored dBG : over {} samples and has {} k-mers and {} color-classes.",
					cdbg.get_num_samples(), cdbg.dbg.dist_elts(), colorClassCount);

	console -> info("Total time taken = {} s.", t_end - t_start);

// TODO what the hell should we do
//	console -> info("Merged CQF metadata:");
//	cdbg.dbg.dump_metadata();


	console->info("Done with cqf merge");
//    uint64_t num_colorBuffers = 4;
//    MSTQuery mst1(cdbg1.prefix, 23, 23, 10, console);
//    MSTQuery mst2(cdbg2.prefix, 23, 23, 10, console);

//    mst1.storeStructure();
//    mst2.storeStructure();
	console->info("{}, {}", cdbg1.prefix, cdbg2.prefix);
	MSTMerger mst(&cdbg.dbg, cdbg.prefix, console, threadCount, cdbg1.prefix, cdbg2.prefix, num_colorBuffers);
	console->info("MST Initiated. Now merging the two MSTs..");
	mst.mergeMSTs();
	serialize_cqf_and_sampleid_list();


}

#endif
