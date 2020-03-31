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

#include "xxhash.h"
#include "BooPHF.h"

#include "mstMerger.h"
//#include <thread>
#include <future>
#include <unistd.h>

// Required to hash colo-id pair objects. Resorted to boost::hash_combine
// instead of plain XOR hashing. For more explanation, consult
// https://stackoverflow.com/questions/35985960/c-why-is-boosthash-combine-the-best-way-to-combine-hash-values
class Custom_Pair_Hasher {
public:
    uint64_t operator()(const std::pair<colorIdType, colorIdType> &key, uint64_t seed = 0) const {
        seed ^= std::hash<uint64_t>{}(key.first) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
        seed ^= std::hash<uint64_t>{}(key.second) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
        return seed;
    }
};

struct ColorPair {
    colorIdType c1, c2;

    ColorPair(colorIdType c1, colorIdType c2) : c1(c1), c2(c2) {}
    ColorPair() {c1=0; c2=0;}

    bool operator==(const ColorPair &rhs) const {
        return c1 == rhs.c1 and c2 == rhs.c2; }

    bool operator!=(const ColorPair &rhs) const {
        return c1 != rhs.c1 or c2 != rhs.c2; }

};

class Fat_Custom_Pair_Hasher {
public:
    uint64_t operator()(const ColorPair&key, uint64_t seed = 0) const {
        return XXH64(reinterpret_cast<void*>(const_cast<ColorPair*>(&key)), sizeof(key), seed);
    }
};

typedef boomphf::mphf<ColorPair, Fat_Custom_Pair_Hasher> boophf_t;

// Hash-map type for color-id pairs.
typedef std::unordered_map<std::pair<colorIdType, colorIdType>, uint64_t, Custom_Pair_Hasher> idPairMap_t;

// adapted from :
// http://stackoverflow.com/questions/34875315/implementation-my-own-list-and-iterator-stl-c
class ColorIdPairIterator {
public:
    using self_type = ColorIdPairIterator;
    using value_type = ColorPair;//std::pair<colorIdType, colorIdType>;
    using reference = value_type &;
    using pointer = value_type *;
    using iterator_category = std::forward_iterator_tag;
    using difference_type = int64_t;

//    ColorIdPairIterator() = delete;

    ColorIdPairIterator(std::string &inputFile, bool isEnd = false) {
        input_.open(inputFile);
        curIt = std::istream_iterator<colorIdType>(input_);
        end_.open(inputFile);
        end_.seekg(0, std::ifstream::end);
        endIt = std::istream_iterator<colorIdType >(end_);
        if (isEnd) {
            curIt = endIt;
        }
        std::cerr << (*curIt) << " vs " << (*endIt) << "\n";
        if (curIt != endIt) {
            for (auto i = 0; i < 20; i++) {
                std::cerr << (*curIt) << "\n";
                curIt++;
            }
            advance_();
        }
    }

    ColorIdPairIterator(const ColorIdPairIterator &other) {
//        fileName = other.fileName;
//        input_.open(fileName);
//        input_.seekg(other.input_.tellg());
        val_ = other.val_;
        curIt = other.curIt;
        endIt = other.endIt;
    }

    ColorIdPairIterator operator++() {
        ColorIdPairIterator i = *this;
        advance_();
        return i;
    }

    const ColorIdPairIterator operator++(int) {
        advance_();
        return *this;
    }

    reference operator*() {
        return val_;
    }

    pointer operator->() {
         return &val_;
    }

    bool operator==(const self_type &rhs) {
        return curIt == rhs.curIt and val_ == rhs.val_; }

    bool operator!=(const self_type &rhs) {
        return curIt != rhs.curIt or val_ != rhs.val_; }

private:
    void advance_() {
        if (curIt == endIt) {
            val_ = ColorPair(static_cast<colorIdType >(invalid), static_cast<colorIdType >(invalid));
            std::cerr << "end: " << val_.c1 << " " << val_.c2 << "\n";
            return;
        }
        auto c1 = *curIt;
        curIt++;
        auto c2 = *curIt;
        val_ = ColorPair(c1, c2);
        curIt++;
    }

//    std::string fileName;
    std::ifstream input_;
    std::ifstream end_;
    ColorPair val_;
    std::istream_iterator<colorIdType> curIt;
    std::istream_iterator<colorIdType> endIt;
};


template<class qf_obj, class key_obj>
class CQF_merger {
public:
    CQF_merger(std::string &firstCQF, std::string &secondCQF, std::string &outputCQF, spdlog::logger *logger,
               uint64_t threadNum);

    // Merges two Colored dBG 'cdbg1' and 'cdbg2' into the colored dBG 'cdbg'
    // (all the CdBG's are private members).
    void merge();


private:

    // Name of the temporary list of color-id pairs.
    const std::string EQ_ID_PAIRS_FILE = std::string("color-id-pairs");

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

    // Blocks for minimizers in the output CDBG
    std::vector<uint64_t> minimizerBlocks;
    std::vector<std::unique_ptr<std::vector<std::pair<uint64_t, colorIdType>>>> minimizerKeyColorList[2];
    uint64_t kbits;
    qf_hashmode hashmode;
    uint32_t seed;
    uint64_t kmerMask;



    // Hash-map for the sampled (on abundance) color-id pairs.
    // Used as the form (pair -> abundance) earlier, and finally as (pair -> colorId).
    idPairMap_t sampledPairs;

    std::unique_ptr<boophf_t> colorMph;

    // Samples 'SAMPLE_PAIR_COUNT' number of most abundant color-id pairs from the
    // first 'sampleKmerCount' distinct k-mers of the CdBGs 'cdbg1' and 'cdbg2',
    // into the map 'sampledPairs', which is of the format (pair -> abundance).
    uint64_t sample_color_id_pairs(uint64_t sampleKmerCount);

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
    uint64_t fill_disk_bucket(uint64_t startingBlock = 0);

    // Filters the disk-buckets to contain only unique color-id pairs.
    void filter_disk_buckets();

    // Builds an MPH (Minimal Perfect Hash) table for each disk-bucket.
    void build_MPH_tables();

    // Builds the output merged CQF.
    void build_CQF();

    // Given a color-id pair 'idPair', returns the newly assigned color-id of this
    // pair at the merged CdBG.
    inline uint64_t get_color_id(const std::pair<uint64_t, uint64_t> &idPair);

    // Serializes the output CQF and sample-id mapping.
    void serializeRemainingStructures();

    void store_color_pairs();

    uint64_t walkBlockedCQF(ColoredDbg<qf_obj, key_obj> &curCdbg, uint64_t curBlock, bool isFirst);
};


#endif
