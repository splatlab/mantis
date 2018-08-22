//
// Created by Fatemeh Almodaresi on 8/17/18.
//

#include "colorEncoder.h"
#include <iostream>
#include <iterator>
#include <algorithm>

bool ColorEncoder::addColorClass(uint64_t kmer, uint64_t eqId, const std::vector<uint32_t> bv) {
    // create list of edges to be processed
    // 1. list of neighbors
    // 2. zero to node
    // calc. the distance between the node and any of the neighbors
    // that exist and the edge hasn't been seen
    dna::canonical_kmer cur(k, kmer);
    std::unordered_set<std::pair<uint64_t, uint64_t>, pair_hash> newEdges;
    // case 1. edges between the node and its neighbors
    for (auto nei_eqId : neighbors(cur)) {
        uint64_t cur_eqId{eqId};
        if (nei_eqId != cur_eqId) {
            if (nei_eqId < cur_eqId) {
                std::swap(nei_eqId, cur_eqId);
            }
            if (!hasEdge(cur_eqId, nei_eqId)) {
                newEdges.insert(std::make_pair(cur_eqId, nei_eqId));
            }
        }
    }
    // case 2. edge from zero to the node
    if (!hasEdge(zero, eqId)) {
        newEdges.insert(std::make_pair(zero, eqId));
    }
    for (auto &newEdge : newEdges) {
        auto deltas = hammingDist(newEdge.first, newEdge.second);
        updateMST(newEdge.first, newEdge.second, deltas);
        addEdge(newEdge.first, newEdge.second, deltas.size());
    }
    if (newEdges.size())
        return true;
    return false;
}

// deltas should *NOT* be passed by reference
bool ColorEncoder::updateMST(uint64_t n1, uint64_t n2, std::vector<uint64_t> deltas) { // n2 > n1
    // If n2 is a new color Id (next color Id)
    if (n2 == colorClsCnt) {
        parentbv[n2] = n1;
        deltaM.insertDeltas(n2, deltas);
        colorClsCnt++;
        if (parentbv.size() == colorClsCnt) {
            parentbv.resize(parentbv.size()+parentbv.size()/2);
        }
        return true;
    }
    // find the max weight edge from each of the ends to their lca, called w1, w2
    uint64_t w1, w2, p1, p2;
    uint64_t w = deltas.size();
    maxWeightsTillLCA(n1, n2, w1, w2, p1, p2);
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
    while (child != p2) {
        auto tmp = parentbv[child];
        parentbv[child] = parent;
        auto prevDeltas = deltaM.getDeltas(child);
        deltaM.insertDeltas(child, deltas);
        deltas = prevDeltas;
        parent = n2;
        child = tmp;
    }
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
    while (i != zero) {
        if (lru_cache.contains(i)) {
            const auto &vs = lru_cache[i];
            for (auto v : vs) {
                xorflips[v] = 1;
            }
            foundCache = true;
            break;
        }
        deltaIndices.push_back(i);
        i = iparent;
        iparent = parentbv[i];
    }

    uint64_t pctr{0};
    for (auto index : deltaIndices) {
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

std::vector<uint64_t> ColorEncoder::hammingDist(uint64_t i, uint64_t j) {
    std::vector<uint64_t> res{0};
    auto n1 = buildColor(i);
    auto n2 = buildColor(j);
    // merge
    // with slight difference of not inserting values that appear in both vectors
    uint64_t i1{0}, i2{0};
    while (i1 < n1.size() or i2 < n2.size()) {
        if (i1 == n1.size()) {
            copy(n2.begin()+i2, n2.end(), back_inserter(res));
            i2 = n2.size();
        }
        else if (i2 == n2.size()) {
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
    return res;
}

void ColorEncoder::maxWeightsTillLCA(uint64_t n1, uint64_t n2,
                                     uint64_t &w1, uint64_t &w2,
                                     uint64_t &p1, uint64_t &p2) {
}

std::unordered_set<uint64_t> ColorEncoder::neighbors(dna::canonical_kmer n) {
    std::unordered_set<uint64_t > result;
    for (const auto b : dna::bases) {
        uint64_t eqid, idx;
        if (exists(b >> n, eqid))
            result.insert(eqid);
        if (exists(n << b, eqid))
            result.insert(eqid);
    }
    return result;
}

bool ColorEncoder::exists(dna::canonical_kmer e, uint64_t &eqid) {
    uint64_t tmp = e.val;
    KeyObject key(HashUtil::hash_64(tmp, BITMASK(cqf.keybits())), 0, 0);
    auto eq_idx = cqf.query(key);
    return eq_idx != 0;
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

int main(int argc, char* argv[]) {

}