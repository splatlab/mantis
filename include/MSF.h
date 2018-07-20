//
// Created by Fatemeh Almodaresi on 7/20/18.
//

#ifndef MANTIS_MSF_H
#define MANTIS_MSF_H
#include<bits/stdc++.h>
#include <sstream>
#include <unordered_set>
#include <queue>
#include "clipp.h"
#include "bitvector.h"
//#include "sdsl/bits.hpp"

#define EQS_PER_SLOT 20000000

using namespace std;

typedef std::vector<sdsl::rrr_vector < 63>> eqvec;

struct Edge {
    uint32_t n1;
    uint32_t n2;
    uint16_t weight;

    Edge(uint32_t inN1, uint32_t inN2, uint16_t inWeight)
            : n1(inN1), n2(inN2), weight(inWeight) {}
};

struct EdgePtr {
    uint16_t bucket;
    uint32_t idx;

    EdgePtr(uint16_t bucketIn, uint32_t idxIn) : bucket(bucketIn), idx(idxIn) {}
};

struct Child {
    uint32_t id;
    uint16_t weight;

    Child(uint32_t inN1, uint16_t inWeight) : id(inN1), weight(inWeight) {}
};

struct Path {
    uint32_t id;
    uint32_t steps;
    uint64_t weight;

    Path(uint32_t idIn,
         uint32_t stepsIn,
         uint64_t weightIn) : id(idIn), steps(stepsIn), weight(weightIn) {}
};

struct DisjointSetNode {
    uint32_t parent{0};
    uint64_t rnk{0}, w{0}, edges{0};

    void setParent(uint32_t p) { parent = p; }

    void mergeWith(DisjointSetNode &n, uint16_t edgeW, uint32_t id) {
        n.setParent(parent);
        w += (n.w + static_cast<uint64_t>(edgeW));
        edges += (n.edges + 1);
        n.edges = 0;
        n.w = 0;
        if (rnk == n.rnk) {
            rnk++;
        }
    }
};

// To represent Disjoint Sets
struct DisjointSets {
    std::vector<DisjointSetNode> els;
    uint64_t n;

    // Constructor.
    DisjointSets(uint64_t n) {
        // Allocate memory
        this->n = n;
        els.resize(n);
        // Initially, all vertices are in
        // different sets and have rank 0.
        for (uint64_t i = 0; i <= n; i++) {
            //every element is parent of itself
            els[i].setParent(i);
        }
    }

    // Find the parent of a node 'u'
    // Path Compression
    uint32_t find(uint32_t u) {
        /* Make the parent of the nodes in the path
           from u--> parent[u] point to parent[u] */
        if (u != els[u].parent)
            els[u].parent = find(els[u].parent);
        return els[u].parent;
    }

    // Union by rank
    void merge(uint32_t x, uint32_t y, uint16_t edgeW) {
        x = find(x), y = find(y);

        /* Make tree with smaller height
           a subtree of the other tree  */
        if (els[x].rnk > els[y].rnk) {
            els[x].mergeWith(els[y], edgeW, x);

        } else {// If rnk[x] <= rnk[y]
            els[y].mergeWith(els[x], edgeW, y);
        }
    }
};

// Structure to represent a graph
struct Graph {

    uint64_t V;

    vector<vector<Edge>> edges;
    vector<vector<EdgePtr>> mst;

    uint64_t mst_totalWeight{0};

    Graph(uint64_t bucketCnt) { edges.resize(bucketCnt); }

    // Utility function to add an edge
    void addEdge(uint32_t u, uint32_t v, uint16_t w) {
        edges[w - 1].emplace_back(u, v, w);
        //edges.emplace_back(u, v, w);
    }

    // Function to find MST using Kruskal's
    // MST algorithm
    DisjointSets kruskalMSF(uint32_t bucketCnt) {
        int mst_wt = 0; // Initialize result

        // Create disjoint sets
        DisjointSets ds(V);

        std::string tmp;
        uint64_t n1{0}, n2{0}, cntr{0}, mergeCntr{0};
        uint32_t w{0};
        sdsl::bit_vector nodes(V, 0);
        // Iterate through all sorted edges
        for (auto bucketCntr = 0; bucketCntr < bucketCnt; bucketCntr++) {
            //ifstream file(filename);
            /*std::getline(file, tmp);
            while (file.good()) {
                file >> n1 >> n2 >> w;*/
            uint32_t edgeIdxInBucket = 0;
            for (auto it = edges[bucketCntr].begin(); it != edges[bucketCntr].end(); it++) {
                //if (w == bucketCntr) {
                w = it->weight;
                uint32_t u = it->n1;
                uint32_t v = it->n2;
                uint32_t set_u = ds.find(u);
                uint32_t set_v = ds.find(v);

                // Check if the selected edge is creating
                // a cycle or not (Cycle is created if u
                // and v belong to same set)
                if (set_u != set_v) {
                    // Current edge will be in the MST
                    // Merge two sets
                    ds.merge(set_u, set_v, w);
                    mst[u].emplace_back(bucketCntr, edgeIdxInBucket);
                    mst[v].emplace_back(bucketCntr, edgeIdxInBucket);
                    nodes[u] = 1;
                    nodes[v] = 1;
                    mst_totalWeight += w;
                    mergeCntr++;
                }/* else {
                            if (nodes.find(u) == nodes.end() || nodes.find(v) == nodes.end())
                                std::cerr << u << " " << v << " " << set_u << " " << set_v << "\n";
                    }*/
                cntr++;
                if (cntr % 1000000 == 0) {
                    std::cerr << "edge " << cntr << " " << mergeCntr << "\n";
                }
                edgeIdxInBucket++;
                //}
            }
            /*file.clear();
            file.seekg(0, file.beg);*/

        }
        //file.close();
        uint64_t distinctNodes{0};
        for (uint64_t i = 0; i < V; i += 64) {
            distinctNodes += sdsl::bits::cnt(nodes.get_int(i, 64));
        }

        std::cerr << "final # of edges: " << cntr
                  << "\n# of merges: " << mergeCntr
                  << "\n# of distinct nodes: " << distinctNodes
                  << "\n";
        return ds;
    }
};

void loadEqs(std::string filename, eqvec &bvs) {
    bvs.reserve(20);
    std::string eqfile;
    std::ifstream eqlist(filename);
    if (eqlist.is_open()) {
        uint64_t accumTotalEqCls = 0;
        while (getline(eqlist, eqfile)) {
            sdsl::rrr_vector<63> bv;
            bvs.push_back(bv);
            sdsl::load_from_file(bvs.back(), eqfile);
        }
    }
    std::cerr << "loaded all the equivalence classes: "
              << ((bvs.size() - 1) * EQS_PER_SLOT + bvs.back().size())
              << "\n";
}

void buildColor(eqvec &bvs,
                std::vector<uint64_t> &eq,
                uint64_t eqid,
uint64_t num_samples) {
uint64_t i{0}, bitcnt{0}, wrdcnt{0};
uint64_t idx = eqid / EQS_PER_SLOT;
uint64_t offset = eqid % EQS_PER_SLOT;
//std::cerr << eqid << " " << num_samples << " " << idx << " " << offset << "\n";
while (i<num_samples) {
bitcnt = std::min(num_samples - i, (uint64_t) 64);
uint64_t wrd = (bvs[idx]).get_int(offset * num_samples + i, bitcnt);
eq[wrdcnt++] = wrd;
i += bitcnt;
}
}

uint16_t sum1s(eqvec &bvs, uint64_t eqid,
uint64_t num_samples, uint64_t numWrds) {
uint16_t res{0};
std::vector<uint64_t> eq;
eq.resize(numWrds);
buildColor(bvs, eq, eqid, num_samples);
for (uint64_t i = 0; i < eq.size(); i += 1) {
res += (uint16_t)sdsl::bits::cnt(eq[i]);
}
return res;
}

// for two non-zero nodes, delta list is positions that xor of the bits was 1
std::vector<uint32_t> getDeltaList(eqvec &bvs,
                                   uint64_t eqid1,uint64_t eqid2, uint64_t num_samples, uint64_t numWrds) {
std::vector<uint32_t> res;
std::vector<uint64_t> eq1, eq2;
eq1.resize(numWrds);
eq2.resize(numWrds);
buildColor(bvs, eq1, eqid1, num_samples);
buildColor(bvs, eq2, eqid2, num_samples);

for (uint32_t i = 0; i < eq1.size(); i += 1) {
uint64_t eq12xor = eq1[i] ^ eq2[i];
for (uint32_t j = 0; j < 64; j++) {
if ( (eq12xor >> j) & 0x01 ) {
res.push_back(i*64+j);
}
}
}

return res; // rely on c++ optimization
}

// for those connected to node zero, delta list is position of set bits
std::vector<uint32_t> getDeltaList(eqvec &bvs,
                                   uint64_t eqid1, uint64_t num_samples, uint64_t numWrds) {
std::vector<uint32_t> res;
std::vector<uint64_t> eq1;
eq1.resize(numWrds);
buildColor(bvs, eq1, eqid1, num_samples);

for (uint32_t i = 0; i < eq1.size(); i += 1) {
for (uint32_t j = 0; j < 64; j++) {
if ( (eq1[i] >> j) & 0x01 ) {
res.push_back(i*64+j);
}
}
}

return res; // rely on c++ optimization
}

#endif //MANTIS_MSF_H
