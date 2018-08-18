//
// Created by Fatemeh Almodaresi on 8/17/18.
//

#include "ColorEncoder.h"
#include <iostream>

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
    DeltaManager deltaManager(3000, 10, 6);
    std::vector<uint64_t> d;
    d.push_back(0);
    d.push_back(300);
    d.push_back(700);
    deltaManager.insertDeltas(0, d);

    d.clear();
    d.push_back(2999);
    deltaManager.insertDeltas(1, d);

    d.clear();
    d.push_back(100);
    d.push_back(200);
    d.push_back(300);
    d.push_back(400);
    d.push_back(500);
    d.push_back(600);
    /*d.push_back(700);
    d.push_back(800);
    d.push_back(900);
    d.push_back(1000);
    d.push_back(1100);*/
    deltaManager.insertDeltas(2, d);

    auto vec = deltaManager.getDeltas(2);
    std::cerr << "\n\ndeltas for 2\n";
    for (auto v : vec) {
        std::cerr << v << " ";
    }

    std::cerr << "\n";
    vec = deltaManager.getDeltas(0);
    std::cerr << "\n\ndeltas for 0\n";
    for (auto v : vec) {
        std::cerr << v << " ";
    }

    vec = deltaManager.getDeltas(1);
    std::cerr << "\n\ndeltas for 1\n";
    for (auto v : vec) {
        std::cerr << v << " ";
    }
/*
    deltaManager.swapDeltas(0,2);
    std::cerr << "After swap:\n";
    vec = deltaManager.getDeltas(0);
    std::cerr << "\n\ndeltas for 0\n";
    for (auto v : vec) {
        std::cerr << v << " ";
    }

    vec = deltaManager.getDeltas(2);
    std::cerr << "\n\ndeltas for 2\n";
    for (auto v : vec) {
        std::cerr << v << " ";
    }
    std::cerr << "\n";*/
}