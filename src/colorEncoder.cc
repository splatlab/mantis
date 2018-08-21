//
// Created by Fatemeh Almodaresi on 8/17/18.
//

#include "colorEncoder.h"
#include <iostream>

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
    if (n2 == colorClsCnt + 1) {
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

std::vector<uint64_t> ColorEncoder::buildColor(uint64_t eqid) {
    std::vector<uint32_t> flips(numSamples);
    std::vector<uint32_t> xorflips(numSamples, 0);
/*
    uint64_t i{eqid}, from{0}, to{0};
    bool foundCache = false;
    uint32_t iparent = parentbv[i];
    while (iparent != i) {
        if (lru_cache.contains(i)) {
            const auto &vs = lru_cache[i];
            for (auto v : vs) {
                xorflips[v] = 1;
            }
            foundCache = true;
            break;
        }
        from = (i > 0) ? (sbbv(i) + 1) : 0;
        froms.push_back(from);

        i = iparent;
        iparent = parentbv[i];
    }
    if (!foundCache and i != zero) {
        from = (i > 0) ? (sbbv(i) + 1) : 0;
        froms.push_back(from);
    }

    uint64_t pctr{0};
    for (auto f : froms) {
        bool found = false;
        uint64_t wrd{0};
        //auto j = f;
        uint64_t offset{0};
        auto start = f;
        do {
            wrd = bbv.get_int(start, 64);
            for (uint64_t j = 0; j < 64; j++) {
                flips[deltabv[start + j]] ^= 0x01;
                if ((wrd >> j) & 0x01) {
                    found = true;
                    break;
                }
            }
            start += 64;
        } while (!found);
    }
*/

    // return the indices of set bits
    std::vector<uint64_t> eq;
    uint64_t numWrds = numSamples/64+1;
    eq.reserve(numWrds);
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
    // If i is the dummy node zero
    // ow
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