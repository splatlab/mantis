//
// Created by Fatemeh Almodaresi on 2019-05-10.
//

#ifndef MANTIS_MHUTIGBUILDER_H
#define MANTIS_MHUTIGBUILDER_H

#include <set>

#include "mantisconfig.hpp"
#include "kmer.h"
#include "gqf/hashutil.h"
#include "gqf_cpp.h"
#include "spdlog/spdlog.h"
#include "canonicalKmer.h"

typedef uint32_t colorIdType;

class MHUtigBuilder {
public:
    MHUtigBuilder(std::string prefixIn, std::shared_ptr<spdlog::logger> loggerIn, uint32_t numThreads);
    std::set<dna::base> prevNeighbors(CQF<KeyObject> &cqf, dna::kmer& n);
    std::set<dna::base> nextNeighbors(CQF<KeyObject> &cqf, dna::kmer& n);
    bool exists(CQF<KeyObject> &cqf, dna::kmer e, colorIdType &eqid);
    void buildUnitigs(CQF<KeyObject> &cqf);
private:
    uint32_t k;
    spdlog::logger* logger;
    uint32_t nThreads;
    std::string prefix;
    CQF<KeyObject>* cqf;

};
#endif //MANTIS_MHUTIGBUILDER_H
