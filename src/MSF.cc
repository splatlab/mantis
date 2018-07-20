//
// Created by Fatemeh Almodaresi on 7/11/18.
//
// The algorithm's basic implementation taken from
// https://www.geeksforgeeks.org/kruskals-minimum-spanning-tree-using-stl-in-c/
//

#include<bits/stdc++.h>
#include <sstream>
#include <unordered_set>
#include <queue>
#include "clipp.h"
#include "bitvector.h"
//#include "sdsl/bits.hpp"

#define EQS_PER_SLOT 20000000

using namespace std;

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

struct Opts {
    std::string filename;
    uint64_t numNodes; // = std::stoull(argv[2]);
    uint32_t bucketCnt; // = std::stoull(argv[3]);
    uint16_t numSamples;
    std::string eqClsListFile;
};

void loadEqs(std::string filename, std::vector<sdsl::rrr_vector < 63>>

&bvs) {
bvs.reserve(20);
std::string eqfile;
std::ifstream eqlist(filename);
if (eqlist.

is_open()

) {
uint64_t accumTotalEqCls = 0;
while (
getline(eqlist, eqfile
)) {
sdsl::rrr_vector<63> bv;
bvs.
push_back(bv);
sdsl::load_from_file(bvs
.

back(), eqfile

);
}
}
//BitVectorRRR bv(eqfile);
std::cerr << "loaded all the equivalence classes: "
<< ((bvs.

size()

- 1) * EQS_PER_SLOT + bvs.

back()

.

size()

)
<< "\n";
}

void buildColor(std::vector<sdsl::rrr_vector < 63>> &bvs,
std::vector<uint64_t> &eq,
        uint64_t eqid,
uint64_t num_samples
) {
uint64_t i{0}, bitcnt{0}, wrdcnt{0};
uint64_t idx = eqid / EQS_PER_SLOT;
uint64_t offset = eqid % EQS_PER_SLOT;
//std::cerr << eqid << " " << num_samples << " " << idx << " " << offset << "\n";
while (i<num_samples) {
bitcnt = std::min(num_samples - i, (uint64_t) 64);
uint64_t wrd = (bvs[idx]).get_int(offset * num_samples + i, bitcnt);
eq[wrdcnt++] =
wrd;
i +=
bitcnt;
}
}

uint16_t sum1s(std::vector<sdsl::rrr_vector < 63>> &bvs, uint64_t eqid,
        uint64_t num_samples) {
    uint16_t res{0};
    std::vector<uint64_t> eq;
    eq.resize((uint64_t)std::ceil((double)num_samples/64.0));
    buildColor(bvs, eq, eqid, num_samples);
    for (uint64_t i = 0; i < eq.size(); i += 1) {
        res += (uint16_t)sdsl::bits::cnt(eq[i]);
    }
    return res;
}

int main(int argc, char *argv[]) {
    /* Let us create above shown weighted
       and undirected graph */
    using namespace clipp;
    enum class mode {
        build, help
    };
    mode selected = mode::help;
    Opts opt;

    auto build_mode = (
            command("build").set(selected, mode::build),
                    required("-e", "--edge-filename") &
                    value("edge_filename", opt.filename) % "File containing list of eq. class edges.",
                    required("-n", "--eqCls-cnt") &
                    value("equivalenceClass_count", opt.numNodes) % "Total number of equivalence (color) classes.",
                    required("-b", "--bucket-cnt") &
                    value("bucket_count", opt.bucketCnt) % "Total number of valid distances.",
                    required("-s", "--numSamples") &
                    value("numSamples", opt.numSamples) % "Total number of experiments (samples).",
                    required("-c", "--eqCls-lst") &
                    value("eqCls_list", opt.eqClsListFile) % "File containing list of equivalence (color) classes."
    );

    auto cli = (
            (build_mode | command("help").set(selected, mode::help)
            )
    );

    decltype(parse(argc, argv, cli)) res;
    try {
        res = parse(argc, argv, cli);
    } catch (std::exception &e) {
        std::cout << "\n\nParsing command line failed with exception: " << e.what() << "\n";
        std::cout << "\n\n";
        std::cout << make_man_page(cli, "MSF");
        return 1;
    }
    if (!res) {
        std::cerr << "Cannot parse the input arguments\n";
        std::exit(1);
    }
    if (selected == mode::help) {
        std::cerr << make_man_page(cli, "MSF");
        std::exit(1);
    }
    //explore_options_verbose(res);

    std::cerr << "here are the inputs: \n"
              << opt.filename << "\n"
              << opt.numNodes << "\n"
              << opt.bucketCnt << "\n";
    opt.numNodes++; // number of nodes is one more than input including the zero node

    std::cerr << "Loading all the equivalence classes first .. \n";
    std::vector<sdsl::rrr_vector<63>> eqs;
    eqs.reserve(20);
    loadEqs(opt.eqClsListFile, eqs);
    std::cerr << "Done loading list of equivalence class buckets\n";

    ifstream file(opt.filename);
    Graph g(opt.bucketCnt);

    uint16_t w_;
    uint32_t n1, n2, edgeCntr{0}, zero{(uint32_t) opt.numNodes - 1};
    {
        std::cerr << "Adding edges from 0 to each node with number of set bits in the node as weight .. \n";
        for (uint32_t i = 0; i < opt.numNodes; i++) {
            uint16_t ones = sum1s(eqs, i, opt.numSamples);
            g.addEdge(zero, i, ones);
        }
        std::cerr << "Done adding 0-end edges.\n";
        std::cerr << "Adding edges between color classes .. \n";
        while (file.good()) {
            file >> n1 >> n2 >> w_;
            g.addEdge(n1, n2, w_);
            edgeCntr++;
        }
        //file.clear();
        //file.seekg(0, ios::beg);
        file.close();
        std::cerr << "Done adding edges between color classes .. \n";

        std::cerr << "\n# of edges: " << edgeCntr
                  << "\n";
//        nodes.clear();
    }
    g.V = opt.numNodes;
    g.mst.resize(opt.numNodes);
    //g.V = numNodes;
    //ifstream file(filename);

    DisjointSets ds = g.kruskalMSF(opt.bucketCnt);
    std::queue<uint32_t> q;
    vector<vector<Child>> p2c(opt.numNodes);
    vector<uint32_t> c2p(opt.numNodes);

    for (auto e = 0; e < g.mst.size(); e++) {
        // put all the mst leaves into the queue
        if (g.mst[e].size() == 1) {
            q.push(e); // e is a leaf
        }
    }
    // now go run the algorithm to find the root of the tree and all the edge directions
    uint32_t root;
    bool check = false;
    uint64_t nodeCntr{0};
    while (!q.empty()) {
        uint32_t node = q.front();
        q.pop();
        if (g.mst[node].size() == 0) {
            // this node is the root
            // and this should be the end of the loop
            // just as a validation, I'll continue and expect the while to end here
            if (check) {
                std::cerr << "A TERRIBLE BUG!!!\n"
                          << "finding a node with no edges should only happen once at the root\n";
                std::exit(1);
            }
            check = true;
            root = node;
            continue;
        }
        // what ever is in q (node) is a leaf and has only one edge left
        // fetch the pointer to that edge (bucket+idx)
        auto buck = g.mst[node][0].bucket;
        auto idx = g.mst[node][0].idx;
        // fetch the two ends of the edge and based on the node id decide which one is the src/dest
        // Destination is the one that represents current node which has been selected as leaf sooner than the other end
        uint64_t dest = g.edges[buck][idx].n1;
        uint64_t src = g.edges[buck][idx].n2;
        if (src == node) { // swap src & dest since src is the leaf
            std::swap(src, dest);
        }
        p2c[src].emplace_back(dest, g.edges[buck][idx].weight);
        // erase the edge from src list of edges
        g.mst[dest].erase(std::remove_if(
                g.mst[dest].begin(), g.mst[dest].end(),
                [buck, idx](const EdgePtr &x) {
                    return x.bucket == buck && x.idx == idx;
                }), g.mst[dest].end());
        // erase the edge from dest list of edges
        g.mst[src].erase(std::remove_if(
                g.mst[src].begin(), g.mst[src].end(),
                [buck, idx](const EdgePtr &x) {
                    return x.bucket == buck && x.idx == idx;
                }), g.mst[src].end());
        // the destination has no edges left
        // but if the src has turned to a leaf (node with one edge) after removal of this last edge
        // add it to the queue
        if (g.mst[src].size() == 1) {
            q.push(src);
        }
        nodeCntr++; // just a counter for the log
        if (nodeCntr % 10000000 == 0) {
            std::cerr << nodeCntr << " nodes processed toward root\n";
        }
    }

    // calculate all the stats!!
    // create the data structures
    std::cerr << "Calculate Stats .. \n";
    std::cerr << "Sum of MST weights: " << g.mst_totalWeight << "\n";
    nodeCntr = 0;
    std::queue<Path> pq;
    pq.push(Path(root, 0, 0));
    double avgDegree{0};
    uint64_t internalNodeCnt{0};
    std::cout << "-1\t" << root << "\t" << 0 << "\t" << 0 << "\t" << 0 << "\n";
    while (!pq.empty()) {
        Path p = pq.front();
        pq.pop();
        avgDegree += p2c[p.id].size();
        if (p2c[p.id].size() != 0) {
            internalNodeCnt++;
        }
        for (auto &c : p2c[p.id]) {
            Path cp(c.id, p.steps + 1, p.weight + c.weight);
            std::cout << p.id << "\t"
                      << cp.id << "\t"
                      << cp.steps << "\t"
                      << c.weight << "\t"
                      << cp.weight << "\n";
            pq.push(cp);
        }
        nodeCntr++;
        if (nodeCntr % 10000000 == 0) {
            std::cerr << nodeCntr << " nodes processed from root\n";
        }
    }

    std::cerr << "Sum of MST weights: " << g.mst_totalWeight << "\n"
              << "internal node count: " << internalNodeCnt << "\n"
              << "total out degree: " << avgDegree << "\n"
              << "average degree: " << avgDegree / internalNodeCnt << "\t" << avgDegree / opt.numNodes << "\n";


    return 0;
}

