//
// Created by Fatemeh Almodaresi on 2018-10-15.
//

#ifndef MANTIS_MSTQUERY_H
#define MANTIS_MSTQUERY_H

#include <sparsepp/spp.h>

#include "spdlog/spdlog.h"
#include "sdsl/bit_vectors.hpp"
#include "mantisconfig.hpp"
#include "lru/lru.hpp"
#include "gqf_cpp.h"
#include "common_types.h"
#include "tsl/hopscotch_map.h"
#include "nonstd/optional.hpp"
#include <mutex>
#include "MantisFS.h"
#include "coloreddbg.h"
#include "BooPHF.h"

//using LRUCacheMap =  LRU::Cache<uint64_t, std::shared_ptr<std::vector<uint64_t>>>;
using LRUCacheMap =  LRU::Cache<uint64_t, std::vector<uint64_t>>;

//using LRUCacheMap = HPHP::ConcurrentScalableCache<uint64_t , std::vector<uint64_t >>;

typedef uint64_t colorIdType;

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
    /*std::vector<uint16_t> heightDist;
    std::vector<uint32_t> weightDist;
    std::vector<uint16_t> noCache_heightDist;
    std::vector<uint32_t> noCache_weightDist;*/

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
//    mantis::EqMap cid2expMap;
//    std::vector<std::vector<mantis::ExperimentID>> cid2exp;
    std::string prefix;
    std::string outfile;
    boomphf::mphf<uint64_t, boomphf::SingleHashFunctor<uint64_t>> colorMph;
    uint64_t threadCount = 1;
    bool MSTIsLoaded = false;
public:
    uint32_t queryK;
    uint32_t indexK;
    sdsl::int_vector<> parentbv;
    sdsl::int_vector<> deltabv;
    sdsl::bit_vector::select_1_type sbbv;
    mantis::QueryResults allQueries;

    static colorIdType getNodeCount(std::string &prefix) {
        static sdsl::int_vector<> parentbv;
        sdsl::load_from_file(parentbv, prefix + mantis::PARENTBV_FILE);
        colorIdType size = parentbv.size();
        parentbv.resize(0);
        return size;
    }

    MSTQuery(std::string prefixIn, uint32_t indexKIn, uint32_t queryKIn,
             uint64_t numSamplesIn, spdlog::logger *loggerIn, bool loadMSTAhead=true) :
            numSamples(numSamplesIn), outfile(prefixIn + "queryTmps"), indexK(indexKIn), queryK(queryKIn), logger(loggerIn), prefix(prefixIn) {
        numWrds = (uint64_t) std::ceil((double) numSamples / 64.0);
        if (loadMSTAhead) {
            loadMST(prefix);
            MSTIsLoaded = true;
        }
    }

    MSTQuery(std::string prefixIn, std::string outIn, uint32_t indexKIn, uint32_t queryKIn,
            uint64_t numSamplesIn, spdlog::logger *loggerIn, bool loadMSTAhead=true) :
    numSamples(numSamplesIn), outfile(outIn), indexK(indexKIn), queryK(queryKIn), logger(loggerIn), prefix(prefixIn) {
        numWrds = (uint64_t) std::ceil((double) numSamples / 64.0);
        if (loadMSTAhead) {
            loadMST(prefix);
            MSTIsLoaded = true;
        }
    }

    void loadMST(std::string indexDir);
    std::vector<uint64_t> buildColor(uint64_t eqid, QueryStats &queryStats,
                                     LRUCacheMap *lru_cache,
                                     RankScores* rs,
                                     std::unordered_map<uint64_t, std::vector<uint64_t>> *fixed_cache,
                                     nonstd::optional<uint64_t>& toDecode // output param.  Also decode these
                                     );

    void parseKmers(uint32_t readId, std::string read, uint64_t kmer_size);
    uint64_t parseBulkKmers(std::string &file, u_int64_t kmer_size);
    void findSamples(ColoredDbg<SampleObject<CQF<KeyObject> *>, KeyObject> &cdbg,
                                        LRUCacheMap &lru_cache,
                                        RankScores *rs,
                                        QueryStats &queryStats,
                                        uint64_t numQueries);
//    mantis::QueryResults getResultList(uint64_t numQueries);

    void reset();

    void clear();
    void storeStructure() {
        auto rootID = parentbv.size()-1;
        //BFS
        uint64_t cntr=0;
        std::map<uint64_t, uint64_t> depthMap;
        depthMap[rootID] = 0;
        while (cntr != parentbv.size()-1) {
            uint64_t depth = 0;
            if (depthMap.find(parentbv[cntr]) == depthMap.end()) {
                auto id = cntr;
                auto parentId = parentbv[id];
                while (id != parentId) {
                    depth++;
                    id = parentbv[cntr];
                    parentId = parentbv[id];
                }
            } else {
                depth = depthMap[parentbv[cntr]] + 1;
            }
            depthMap[cntr] = depth;
            cntr++;
            if (cntr % 100000 == 0)
                std::cerr << "\r" << cntr;
        }
        std::cerr << "\n";
        std::ofstream st(prefix + "/mstStructure.out");
        for (auto &kv : depthMap) {
            st << kv.first << " " << kv.second << "\n";
        }
        st.close();
    }

    uint64_t getNumOfDistinctKmers() {
        return kmer2cidMap.size();
    }

    /*void setFixed_size(uint64_t fixed_sizeIn) {
        fixed_size = fixed_sizeIn;
    }*/
};

#endif //MANTIS_MSTQUERY_H
