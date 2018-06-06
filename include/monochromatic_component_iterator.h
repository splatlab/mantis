//
// Created by Fatemeh Almodaresi on 6/4/18.
//

#ifndef MANTIS_MONOCHROME_ITERATOR_H
#define MANTIS_MONOCHROME_ITERATOR_H

#include<set>
#include<queue>
#include <string>
#include <map>

#include <iostream>
#include <fstream>
#include <unordered_map>
#include <set>
#include <unordered_set>

#include <inttypes.h>

#include "sparsepp/spp.h"
#include "tsl/sparse_map.h"
#include "sdsl/bit_vectors.hpp"
//#include "bitvector.h"
#include "kmer.h"
#include "cqf.h"
#include "hashutil.h"
#include "common_types.h"

struct Mc_stats {
    uint64_t nodeCnt = 1;
    uint64_t min_dist = -1; // infinity
};

typedef dna::canonical_kmer edge; // k-mer
typedef dna::canonical_kmer node; // (k-1)-mer

class monochromatic_component_iterator {
public:
    class work_item {
    public:
        node curr;
        uint64_t idx;
        uint64_t colorid;

        work_item(node currin, uint64_t idxin, uint64_t coloridin) : curr(currin), idx(idxin), colorid(coloridin) {}

        bool operator<(const work_item &item2) const {
            return (*this).idx < item2.idx;
        }
    };

    bool done();
    void operator++(void);

    Mc_stats operator*(void);

    monochromatic_component_iterator(const CQF<KeyObject> *g);

private:

    uint32_t k;
    std::queue<work_item> work;


    const CQF<KeyObject> *cqf;
    CQF<KeyObject>::Iterator it;
    sdsl::bit_vector visited;

    bool exists(edge e, uint64_t &idx, uint64_t &eqid);

    std::set<monochromatic_component_iterator::work_item> neighbors(monochromatic_component_iterator::work_item n);

    work_item front(std::queue<work_item> &w);


    uint64_t manhattanDist(uint64_t eq1, uint64_t eq2);

};


#endif //MANTIS_MONOCHROME_ITERATOR_H
