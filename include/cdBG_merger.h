/*
 * ============================================================================
 *
 *         Authors:
 *         Fatemeh Amodaresi, falmodar@cs.umd.edu
 *         Jamshed Khan, jamshed@cs.umd.edu
 *   Organization:  University of Maryland
 *
 * ============================================================================
 */

#ifndef _CDBG_MERGER_H_
#define _CDBG_MERGER_H_

#include "spdlog/spdlog.h"
#include "MantisFS.h"
#include "coloreddbg.h"
#include "BooPHF.h"

#include "mstMerger.h"
//#include <thread>
#include <future>
#include <unistd.h>

template<class qf_obj, class key_obj>
class CdBG_merger {
public:
    CdBG_merger(ColoredDbg<qf_obj, key_obj> &&cdbg1, ColoredDbg<qf_obj, key_obj> &&cdbg2,
                ColoredDbg<qf_obj, key_obj> &&cdbgOut);

    void set_console(spdlog::logger *c) { console = c; }

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
//		const static uint64_t SAMPLE_PAIR_COUNT = std::min((uint64_t)1000000, mantis::NUM_BV_BUFFER);

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
    spdlog::logger *console;

    // Initial time-stamp of object creation.
    std::time_t start_time_;

    // Number of processor-threads to be used at the intermediate steps of  unique
    // color-id pairs filtering and MPH's building.
    uint threadCount = 1;

    // Utility information to display at the result summary.
    uint64_t colorCount1 = 0, colorCount2 = 0;

    // Blocks for minimizers in the output CDBG
    std::vector<uint64_t> minimizerBlocks;
    std::vector<std::unique_ptr<std::vector<std::pair<uint64_t, colorIdType>>>> minimizerKeyColorList[2];
//		std::vector<std::unordered_map<uint64_t, std::pair<colorIdType , colorIdType>>> minimizerKeyColorList;
    uint64_t kbits;
    qf_hashmode hashmode;
    uint32_t seed;
    uint64_t kmerMask;


    // Required to hash colo-id pair objects. Resorted to boost::hash_combine
    // instead of plain XOR hashing. For more explanation, consult
    // https://stackoverflow.com/questions/35985960/c-why-is-boosthash-combine-the-best-way-to-combine-hash-values
    class Custom_Pair_Hasher {
    public:
        uint64_t operator()(const std::pair<uint64_t, uint64_t> &key, uint64_t seed = 0) const {
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

    uint64_t numCCPerBuffer{0};
    uint64_t numCCPerBuffer1{0};
    uint64_t numCCPerBuffer2{0};

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
    inline void add_color_id_pair(uint64_t colorID1, uint64_t colorID2,
                                  std::vector<std::vector<std::ofstream>> &diskBucket);

    // Filters the disk-buckets to contain only unique color-id pairs.
    // Returns the count of unique color-id pairs.
    uint64_t filter_disk_buckets();

    // Builds an MPH (Minimal Perfect Hash) table for each disk-bucket.
    void build_MPH_tables();

    // Builds the output merged CQF.
    void build_CQF();

    // Given a color-id pair 'idPair', returns the newly assigned color-id of this
    // pair at the merged CdBG.
    inline uint64_t get_color_id(const std::pair<uint64_t, uint64_t> &idPair);

    // Serializes the output CQF and sample-id mapping.
    void serializeRemainingStructures();

    // Builds the output color-class bitvectors for the color-id pairs.
    void store_color_pairs(ColoredDbg<qf_obj, key_obj> &cdbg1, ColoredDbg<qf_obj, key_obj> &cdbg2,
                           uint64_t &numColorBuffers);

    uint64_t store_abundant_color_pairs(std::ofstream &output);

    void calc_mst_stats(ColoredDbg<qf_obj, key_obj> &cdbg1, ColoredDbg<qf_obj, key_obj> &cdbg2, std::string &dir1,
                        std::string &dir2);

    uint64_t walkBlockedCQF(ColoredDbg<qf_obj, key_obj> &curCdbg, uint64_t curBlock, bool isFirst);
};


#endif
