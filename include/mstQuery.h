//
// Created by Fatemeh Almodaresi on 2018-10-15.
//

#ifndef MANTIS_MSTQUERY_H
#define MANTIS_MSTQUERY_H

#include "spdlog/spdlog.h"
#include "sdsl/bit_vectors.hpp"
#include "mantisconfig.hpp"
#include "lru/lru.hpp"
#include "gqf_cpp.h"
#include "common_types.h"
#include "tsl/hopscotch_map.h"
#include "nonstd/optional.hpp"

using LRUCacheMap =  LRU::Cache<uint64_t, std::vector<uint64_t>>;

struct QueryStats {
    uint32_t cnt = 0, cacheCntr = 0, noCacheCntr{0};
    uint64_t totSel{0};
    std::chrono::duration<double> selectTime{0};
    std::chrono::duration<double> flipTime{0};
    uint64_t totEqcls{0};
    uint64_t rootedNonZero{0};
    uint64_t nextCacheUpdate{10000};
    uint64_t globalQueryNum{0};
    std::vector<uint64_t> buffer;
    uint64_t numSamples{0};
    tsl::hopscotch_map<uint32_t, uint64_t> numOcc;
    bool trySample{false};
    //std::unordered_map<uint32_t, uint64_t> numOcc;
};

class RankScores {
public:
    RankScores(uint64_t nranks) {rs_.resize(nranks);}

    std::unordered_map<uint32_t, uint32_t>& operator[](uint32_t r) {
        if (r > maxRank_) {
            maxRank_ = std::min(r, static_cast<uint32_t>(rs_.size()-1));
        }
        return (r < rs_.size()) ? rs_[r] : rs_.back();
    }

    void clear() {
        for (auto& m : rs_){ m.clear(); }
        maxRank_ = 0;
    }

    uint32_t maxRank() const { return maxRank_; }

private:
    std::vector<std::unordered_map<uint32_t, uint32_t>> rs_;
    uint32_t maxRank_{0};
};

class MSTQuery {
private:
    uint64_t numSamples;
    uint64_t numWrds;
    uint32_t zero;
    sdsl::bit_vector bbv;
    spdlog::logger *logger{nullptr};
    mantis::QueryMap kmer2cidMap;
    mantis::EqMap cid2expMap;

public:
    uint32_t queryK;
    uint32_t indexK;
    sdsl::int_vector<> parentbv;
    sdsl::int_vector<> deltabv;
    sdsl::bit_vector::select_1_type sbbv;

    MSTQuery(std::string prefix, uint32_t indexKIn, uint32_t queryKIn,
            uint64_t numSamplesIn, spdlog::logger *loggerIn) :
    numSamples(numSamplesIn), indexK(indexKIn), queryK(queryKIn), logger(loggerIn) {
        numWrds = (uint64_t) std::ceil((double) numSamples / 64.0);
        loadIdx(prefix);
    }

    void loadIdx(std::string indexDir);
    std::vector<uint64_t> buildColor(uint64_t eqid, QueryStats &queryStats,
                                     LRUCacheMap *lru_cache,
                                     RankScores* rs,
                                     nonstd::optional<uint64_t>& toDecode // output param.  Also decode these
                                     );

    void parseKmers(std::string read, uint64_t kmer_size);
    void findSamples(CQF<KeyObject> &dbg,
                                        LRUCacheMap &lru_cache,
                                        RankScores *rs,
                                        QueryStats &queryStats);
    mantis::QueryResult convertIndexK2QueryK(std::string &read);

    mantis::QueryResult getResultList();

    void reset();

    uint64_t getNumOfDistinctKmers() {
        return kmer2cidMap.size();
    }
};

#endif //MANTIS_MSTQUERY_H
