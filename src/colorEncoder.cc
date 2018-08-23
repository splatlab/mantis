//
// Created by Fatemeh Almodaresi on 8/17/18.
//

#include "colorEncoder.h"
#include <iostream>
#include <iterator>
#include <algorithm>

bool ColorEncoder::addColorClass(uint64_t kmer, uint64_t eqId, const sdsl::bit_vector &bv) {
    // create list of edges to be processed
    // 1. zero to node
    // 2. list of neighbors
    // calc. the distance between the node and any of the neighbors
    // that exist and the edge hasn't been seen
    duplicated_dna::canonical_kmer cur(k, HashUtil::hash_64i(kmer, BITMASK(cqf.keybits())));
    std::unordered_set<std::pair<uint64_t, uint64_t>, pair_hash> newEdges;

    std::vector<uint64_t> setBits;
    if (!lru_cache.contains(eqId)) {
        setBits = buildColor(bv);
        lru_cache.emplace(eqId, setBits);
    }
    // case 1. edge from zero to the node
    if (!hasEdge(zero, eqId)) {
        //std::cerr << "case1: " << eqId << " " << colorClsCnt << "\n";
        updateMST(zero, eqId, setBits);
        addEdge(zero, eqId, setBits.size());
    }

    // case 2. edges between the node and its neighbors
    for (auto nei_eqId : neighbors(cur)) {
        uint64_t cur_eqId{eqId};
        if (nei_eqId != cur_eqId) {
            if (nei_eqId < cur_eqId) {
                std::swap(cur_eqId, nei_eqId);
            }
            if (!hasEdge(cur_eqId, nei_eqId)) {
                newEdges.insert(std::make_pair(cur_eqId, nei_eqId));
            }
        }
    }
    //std::cerr << "case2 " << eqId << " " << colorClsCnt << " " << newEdges.size() << "\n";
    for (auto &newEdge : newEdges) {
        auto deltas = hammingDist(newEdge.first, newEdge.second);
        if (updateMST(newEdge.first, newEdge.second, deltas)) {
            //std::cerr << "case2 " << eqId << " " << colorClsCnt << "\n";
        }
        addEdge(newEdge.first, newEdge.second, deltas.size());
    }
    if (newEdges.size())
        return true;
    return false;
}

bool ColorEncoder::serialize(std::string prefix) {
    std::string parentbv_file = prefix + "/parent.bv";
    parentbv.resize(colorClsCnt);
    bool parentSuccessfullyStored = sdsl::store_to_file(parentbv, parentbv_file);
    parentbv.resize(0);

    return deltaM.serialize(prefix) and parentSuccessfullyStored;
}

// deltas should *NOT* be passed by reference
bool ColorEncoder::updateMST(uint64_t n1, uint64_t n2, std::vector<uint64_t> deltas) { // n2 > n1
    //std::cerr << "updateMST: n1= " << n1 << " , n2=" << n2 << " ";
    if (n1 > n2) {
        std::swap(n1, n2);
    }
    while (parentbv.size() < colorClsCnt) {
        parentbv.resize(parentbv.size()+parentbv.size()/2);
    }
    // The only time that we will see the edge zero -> n2 is when n2 is observed for the first time
    if (n1 == zero) {
        parentbv[n2] = n1;
        /*for (auto d : deltas) {
            std::cerr << d << " ";
        }
        std::cerr << "\n";
        std::cerr << "insertDeltas 1\n";*/
        deltaM.insertDeltas(n2, deltas);
        if (colorClsCnt < n2)
            colorClsCnt = n2+1; // n2 is the index
        return true;
    }
    // find the max weight edge from each of the ends to their lca, called w1, w2
    uint64_t w1, w2, p1, p2;
    uint64_t w = deltas.size();
    std::pair<Edge, Edge> lr = maxWeightsTillLCA(n1, n2);
    w1 = lr.first.weight;
    w2 = lr.second.weight;
    p1 = lr.first.parent;
    p2 = lr.second.parent;
    // if w >= w1 and w >= w2, happily do nothing
    if (w >= w1 and w >= w2)
        return false;
    if (w1 > w2) {
        std::swap(w1, w2);
        std::swap(p1, p2);
    }
    // Now we know that we're gonna add the edge n1 -> n2
    // and remove edge with weight w2 along the path from n2 to the LCA
    // update parentbv and deltas for n2
    // change the parent/child relationship
    // starting from node n2 toward the c2, child node of the edge with weight w2
    auto parent = n1;
    auto child = n2;
    //6 8 0 0 4
    //std::cerr << "here: " << n1 << " " << n2 << " " << parentbv[n2] << " " << p1 << " " << p2 << " ";
    while (child != p2) {
        uint64_t tmp = parentbv[child];
        //std::cerr << " first: " << tmp << " then ";
        parentbv[child] = parent;
        auto prevDeltas = deltaM.getDeltas(child);
        /*if (n1 == 6 and n2 == 8)
         std::cerr << parent << "->" << child << " " << tmp << ", ";*/
        deltaM.insertDeltas(child, deltas);
        deltas = prevDeltas;
        parent = child;
        child = tmp;
    }
    //std::cerr << "\n";
    return true;
}

// returns list of set bits
std::vector<uint64_t> ColorEncoder::buildColor(uint64_t eqid) {
    std::vector<uint64_t> eq;
    if (eqid == zero) { // if dummy node "zero", return an empty list (no set bits)
        return eq;
    }
    uint64_t numWrds = numSamples/64+1;
    eq.reserve(numWrds);

    std::vector<uint32_t> flips(numSamples);
    std::vector<uint32_t> xorflips(numSamples, 0);
    uint64_t i{eqid};
    std::vector<uint64_t> deltaIndices;
    deltaIndices.reserve(numWrds);
    bool foundCache = false;
    uint32_t iparent = parentbv[i];
    //std::cerr << "\nwalking mst\n";
    while (i != zero) {
        if (lru_cache.contains(i)) {
            const auto &vs = lru_cache[i];
            for (auto v : vs) {
                xorflips[v] = 1;
            }
            foundCache = true;
            break;
        }
        //std::cerr << " " << i;
        deltaIndices.push_back(i);
        i = iparent;
        iparent = parentbv[i];
    }
    //std::cerr << "\n";

    uint64_t pctr{0};
    for (auto index : deltaIndices) {
        //std::cerr << "getDeltas " << index << "\n";
        auto deltas = deltaM.getDeltas(index);
        for (auto d : deltas) {
            flips[d] ^= 0x01;
        }
    }

    // return the indices of set bits
    uint64_t one = 1;
    for (auto i = 0; i < numSamples; i++) {
        if (flips[i] ^ xorflips[i]) {
            eq.push_back(i);
        }
    }
    return eq;
}

std::vector<uint64_t> ColorEncoder::buildColor(const sdsl::bit_vector &bv) {
    std::vector<uint64_t> setBits;
    setBits.reserve(numSamples);
    uint64_t i = 0;
    //std::cerr << " bv bitsize: " << bv.bit_size() << " ";
    while (i < bv.bit_size()) {
        uint64_t bitcnt = numSamples - i >= 64?64:(numSamples - i);
        auto wrd = bv.get_int(i, bitcnt);
        for (uint64_t c=0; c < bitcnt; c++) {
            if ( (wrd >> c) & 0x01) {
                setBits.push_back(i+c);
            }
        }
        i+=64;
    }
    /*std::cerr << " setBits: ";
    for (auto s : setBits) {
        std::cerr << s << " ";
    }
    std::cerr << "\n";*/
    return setBits;
}

std::vector<uint64_t> ColorEncoder::hammingDist(uint64_t i, uint64_t j) {
    //std::cerr << "\nhamming dist between " << i << "," << j << "\n";
    std::vector<uint64_t> res{0};
    auto n1 = buildColor(i);
    auto n2 = buildColor(j);
    //std::cerr << " merge, ";
    // merge
    // with slight difference of not inserting values that appear in both vectors
    uint64_t i1{0}, i2{0};
    while (i1 < n1.size() or i2 < n2.size()) {
        if (i1 == n1.size()) {
            //std::cerr << " i1=n1: " << n1.size() << " " << n2.size() << " ";
            copy(n2.begin()+i2, n2.end(), back_inserter(res));
            i2 = n2.size();
        }
        else if (i2 == n2.size()) {
            //std::cerr << " i2=n2: " << n2.size() << " " << n1.size() << " ";
            copy(n1.begin()+i1, n1.end(), back_inserter(res));
            i1 = n1.size();
        }
        else {
            if (n1[i1] < n2[i2]) {
                res.push_back(n1[i1]);
                i1++;
            } else if (n2[i2] < n1[i1]) {
                res.push_back(n2[i2]);
                i2++;
            } else { //n1[i1] == n2[i2] both have a set bit at the same index
                i1++;
                i2++;
            }
        }
    }
    /*for (auto r : res) {
        std::cerr << r << " ";
    }
    std::cerr << "\n";*/
    return res;
}

//TODO we definitely need to discuss this. This function in its current way is super inefficient
// time/space tradeoff
std::pair<Edge, Edge> ColorEncoder::maxWeightsTillLCA(uint64_t n1, uint64_t n2) {
    std::vector<uint64_t> nodes1;
    std::vector<uint64_t> nodes2;
    uint64_t n = n1;
//    if (n1 == 6 and n2 == 8)
//    std::cerr << "1 ";
    while (n != zero) {
//        if (n1 == 6 and n2 == 8)
//            std::cerr << n << " ";
        nodes1.push_back(n);
        n = parentbv[n];
    }
    nodes1.push_back(zero);
//    if (n1 == 6 and n2 == 8)
//    std::cerr << "\n2 ";
    n = n2;
    while (n != zero) {
//        if (n1 == 6 and n2 == 8)
//            std::cerr << n << " ";
        nodes2.push_back(n);
        n = parentbv[n];
    }
    nodes2.push_back(zero);
//    if (n1 == 6 and n2 == 8)
//    std::cerr << "\n3 lca=";
    auto &n1ref = nodes1;
    auto &n2ref = nodes2;
    if (n1ref.size() < n2ref.size()) {
        std::swap(n1ref, n2ref);
    }

    // find lca
    uint64_t lca = 0;
    for (uint64_t j = 0, i = n1ref.size()-n2ref.size(); j < n2ref.size(); i++, j++) {
        if (n1ref[i] == n2ref[j]) {
            lca = n1ref[i];
//            if (n1 == 6 and n2 == 8)
//            std::cerr << lca;
            break;
        }
    }

    // walk from n1 to lca and find the edge with maximum weight
    Edge e1;
    uint64_t i = 0;
//    if (n1 == 6 and n2 == 8)
//    std::cerr << "\n4 ";
    while (n1ref[i] != lca) {
        auto curW = getEdge(n1ref[i], n1ref[i+1]);
        if (e1.weight < curW) {
//            if (n1 == 6 and n2 == 8)
//            std::cerr << n1ref[i+1] << " " << n1ref[i] << " " << curW << "\n";
            e1 = Edge(n1ref[i+1], n1ref[i], curW);
        }
        i++;
    }

    // walk from n2 to lca and find the edge with maximum weight
    Edge e2;
    i = 0;
    while (n2ref[i] != lca) {
        auto curW = getEdge(n2ref[i], n2ref[i+1]);
        if (e2.weight < curW) {
            e2 = Edge(n2ref[i+1], n2ref[i], curW);
        }
        i++;
    }
    return std::make_pair(e1, e2);
}

std::unordered_set<uint64_t> ColorEncoder::neighbors(duplicated_dna::canonical_kmer n) {
    std::unordered_set<uint64_t > result;
    for (const auto b : duplicated_dna::bases) {
        uint64_t eqid{0}, idx;
        if (exists(b >> n, eqid)) {
            //std::cerr << std::string(n) << "\n" << std::string(b >> n) << "\n";
            result.insert(eqid);
        }
        if (exists(n << b, eqid)) {
            //std::cerr << std::string(n) << "\n" << std::string(b >> n) << "\n";
            result.insert(eqid);
        }
    }
    return result;
}

bool ColorEncoder::exists(duplicated_dna::canonical_kmer e, uint64_t &eqid) {
    uint64_t tmp = e.val;
    KeyObject key(HashUtil::hash_64(tmp, BITMASK(cqf.keybits())), 0, 0);
    eqid = cqf.query(key);
    return eqid != 0;
    // commenting this line, eqIds start from 1 (rather than 0)
/*
        if (eq_idx) {
            eqid = eq_idx - 1; //note be careful about this -1 here. It'll change many things
            return true;
        }
        return false;
*/
}

void ColorEncoder::addEdge(uint64_t i, uint64_t j, uint32_t w) {
    if (i == j) return;
    if (i > j) {
        std::swap(i,j);
    }
    if (edges.find(std::make_pair(i, j)) == edges.end()) {
        edges[std::make_pair(i, j)] = w;
    }
}

bool ColorEncoder::hasEdge(uint64_t i, uint64_t j) {
    if (i > j) {
        std::swap(i, j);
    }
    return i == j or edges.find(std::make_pair(i, j)) != edges.end();
}

uint32_t ColorEncoder::getEdge(uint64_t i, uint64_t j) {
    if (i > j) {
        std::swap(i, j);
    }
    return edges[std::make_pair(i, j)];
}
