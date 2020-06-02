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


struct ColorPair {
    colorIdType c1, c2;

    ColorPair(colorIdType c1, colorIdType c2) : c1(c1), c2(c2) {}
    ColorPair() {c1=0; c2=0;}

    bool operator==(const ColorPair &rhs) const {
        return c1 == rhs.c1 and c2 == rhs.c2; }

    bool operator!=(const ColorPair &rhs) const {
        return c1 != rhs.c1 or c2 != rhs.c2; }

    friend std::istream& operator>> (std::istream& is, ColorPair& cp)
    {
        is>> cp.c1 >> cp.c2;
        return is;
    }

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

struct KeyColorMin {
    uint64_t key;
    colorIdType color;
    uint16_t minimizer;
    KeyColorMin(uint64_t keyIn, colorIdType colorIn, uint16_t minimizerIn): key(keyIn), color(colorIn), minimizer(minimizerIn) {}
};

// adapted from :
// http://stackoverflow.com/questions/34875315/implementation-my-own-list-and-iterator-stl-c
class ColorIdPairIterator {
public:
    using self_type = ColorIdPairIterator;
    using value_type = ColorPair;
    using reference = value_type &;
    using pointer = value_type *;
    using iterator_category = std::forward_iterator_tag;
    using difference_type = int64_t;

//    ColorIdPairIterator() = delete;

    ColorIdPairIterator(std::string &inputFile, bool isEnd = false) {
        fileName = inputFile;
        input_.open(inputFile);
        if (isEnd) {
            input_.seekg(0, std::ios_base::end);
        }
        cntr = 0;
        advance_();
    }

    ColorIdPairIterator(const ColorIdPairIterator &other) {
        fileName = other.fileName;
        input_.open(fileName);
        input_.seekg(other.input_.tellg());
        val_ = other.val_;
        isEnd_ = other.isEnd_;
        cntr = other.cntr;
        buffer = other.buffer;
    }

    ColorIdPairIterator &
    operator=(const ColorIdPairIterator &other) { //}= default;
        fileName = other.fileName;
        input_.open(fileName);
        input_.seekg(other.input_.tellg());
        val_ = other.val_;
        isEnd_ = other.isEnd_;
        cntr = other.cntr;
        buffer = other.buffer;
        return *this;
    }

    const ColorIdPairIterator& operator++() {
        advance_();
        return *this;
    }

  /*  const ColorIdPairIterator& operator++(int) {
        ColorIdPairIterator i = *this;
        advance_();
        return i;
    }*/

    reference operator*() {
        return val_;
    }

    pointer operator->() {
        return &val_;
    }

    bool operator==(const self_type &rhs) {
        return isEnd_ == rhs.isEnd_ and val_ == rhs.val_; }

    bool operator!=(const self_type &rhs) {
        return isEnd_ != rhs.isEnd_ or val_ != rhs.val_; }

    ~ColorIdPairIterator() {
        input_.close();
    }

private:

    void loadIntoBuffer_() {
        auto limit = 1000000;
        buffer.clear();
        buffer.reserve(limit);
        colorIdType c1,c2;
        while (input_.good() and buffer.size() < limit) {
            input_ >> c1;
            if (input_.good()) {
                input_ >> c2;
                buffer.emplace_back(c1, c2);
            }
        }
        cntr = 0;
    }
    void advance_() {
        if (cntr == buffer.size()) {
            loadIntoBuffer_();
        }
        if (not buffer.empty()) {
            val_ = buffer[cntr];
            cntr++;
        } else {
            isEnd_ = true;
            val_ = ColorPair(static_cast<colorIdType >(invalid), static_cast<colorIdType >(invalid));
        }
    }

    std::string fileName;
    mutable std::ifstream input_;
    ColorPair val_;
    bool isEnd_{false};
    uint64_t cntr{0};
    std::vector<ColorPair> buffer;
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
    std::vector<std::unique_ptr<std::vector<uint64_t >>> minimizerKeyList[2];
    std::vector<std::unique_ptr<sdsl::int_vector<>>> minimizerColorList[2];
    uint64_t kbits;
    qf_hashmode hashmode;
    uint32_t seed;
    uint64_t kmerMask;

    uint64_t colorBits[2];
    // Hash-map for the sampled (on abundance) color-id pairs.
    // Used as the form (pair -> abundance) earlier, and finally as (pair -> colorId).
    idPairMap_t sampledPairs;

    std::unique_ptr<boophf_t> colorMph;

    // Samples 'SAMPLE_PAIR_COUNT' number of most abundant color-id pairs from the
    // first 'sampleKmerCount' distinct k-mers of the CdBGs 'cdbg1' and 'cdbg2',
    // into the map 'sampledPairs', which is of the format (pair -> abundance).
    uint64_t sample_colorID_pairs(uint64_t sampleKmerCount);

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
    uint64_t store_colorID_pairs(uint64_t startingBlock = 0);

    // Filters the disk-buckets to contain only unique color-id pairs.
    void sortUniq_colorID_pairs();

    // Builds an MPH (Minimal Perfect Hash) table for each disk-bucket.
    void build_MPHF();

    // Builds the output merged CQF.
    void build_CQF();

    // Given a color-id pair 'idPair', returns the newly assigned color-id of this
    // pair at the merged CdBG.
    inline uint64_t get_colorID(const std::pair<uint64_t, uint64_t> &idPair);

    // Serializes the output CQF and sample-id mapping.
    void serializeRemainingStructures();

    void store_colorID_map();

    uint64_t walkBlockedCQF(ColoredDbg<qf_obj, key_obj> &curCdbg, uint64_t curBlock, bool isFirst);
};


#endif
