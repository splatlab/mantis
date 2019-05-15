//
// Created by Fatemeh Almodaresi on 2019-05-10.
//

#ifndef MANTIS_WALKQFUTIL_H
#define MANTIS_WALKQFUTIL_H

#include <set>

#include "mantisconfig.hpp"
#include "gqf_cpp.h"
#include "spdlog/spdlog.h"

#include "canonicalKmer.h"

typedef uint32_t colorIdType;

struct workItem {
    dna::canonical_kmer node;
    colorIdType colorId;

    workItem(dna::canonical_kmer n, colorIdType c) : node(n), colorId(c) {}

    // Required to be able to use it as a key in set
    bool operator<(const workItem &item2) const {
        return (*this).node < item2.node;
    }
};

class UnitigBuilder {
public:
    UnitigBuilder(std::string prefixIn, std::shared_ptr<spdlog::logger> loggerIn, uint32_t numThreads);
    std::set<workItem> neighbors(CQF<KeyObject> &cqf, workItem n);
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
#endif //MANTIS_WALKQFUTIL_H
