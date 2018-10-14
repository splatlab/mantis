//
// Created by Fatemeh Almodaresi on 2018-10-13.
//

#ifndef MANTIS_MST_H
#define MANTIS_MST_H

#include <set>

// sparsepp include should be above gqf_cpp! ow, we'll get a conflict in MAGIC_NUMBER
#include "sparsepp/spp.h"

#include "mantisconfig.hpp"
#include "gqf_cpp.h"
#include "spdlog/spdlog.h"

#include "canonicalKmer.h"

struct Edge {
    uint32_t n1;
    uint32_t n2;

    Edge(uint32_t inN1, uint32_t inN2) : n1(inN1), n2(inN2) {}

    bool operator==(const Edge& e) const {
        return n1 == e.n1 && n2 == e.n2;
    }
};

// note: @fatal: careful! The hash highly depends on the length of the edge ID (uint32)
struct edge_hash {
    uint64_t operator() (const Edge& e) const {
        uint64_t res = e.n1;
        return (res << 32) | (uint64_t)e.n2;
    }
};

struct workItem {
    dna::canonical_kmer node;
    uint64_t colorId;

    workItem(dna::canonical_kmer n, uint64_t c) : node(n), colorId(c) {}

    // To be able to use in as a key in set
    bool operator<(const workItem &item2) const {
        return (*this).node < item2.node;
    }
};

class MST {
public:
    explicit MST(std::string prefix);
    void buildMST();

private:
    bool buildEdgeSets();
    void findNeighborEdges(CQF<KeyObject>& cqf, CQF<KeyObject>::Iterator it);
    bool calculateWeights();
    bool findMST();

    std::set<workItem> neighbors(CQF<KeyObject>& cqf, workItem n);
    bool exists(CQF<KeyObject>& cqf, dna::canonical_kmer e, uint64_t &eqid);
    inline uint64_t fetchBucketId(uint64_t c1, uint64_t c2);


    std::string prefix;
    uint64_t k;
    uint64_t num_of_ccBuffers;
    uint64_t num_edges = 0;
    std::vector<spp::sparse_hash_set<Edge, edge_hash>> edgesetList;
    std::shared_ptr<spdlog::logger> console{nullptr};
};
#endif //MANTIS_MST_H
