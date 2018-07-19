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
        els.resize(n + 1);
        // Initially, all vertices are in
        // different sets and have rank 0.
        for (uint64_t i = 0; i <= n; i++) {
            //every element is parent of itself
            els[i].setParent(i);
        }
    }

    // Find the parent of a node 'u'
    // Path Compression
    int find(uint32_t u) {
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
    uint64_t V, E;
    //vector<Edge> edges;
    vector<vector<Edge>> edges;
    vector<vector<EdgePtr>> mst;

    Graph(uint64_t bucketCnt) { edges.resize(bucketCnt); }

    // Utility function to add an edge
    void addEdge(uint32_t u, uint32_t v, uint16_t w) {
        edges[w - 1].emplace_back(u, v, w);
        //edges.emplace_back(u, v, w);
    }

    // Function to find MST using Kruskal's
    // MST algorithm
    DisjointSets kruskalMSF(uint32_t bucketCnt/*, ifstream &file*/) {
        int mst_wt = 0; // Initialize result

        // Sort edges in increasing order on basis of cost

        /*sort(edges.begin(), edges.end(),
             [](const Edge &e1, const Edge &e2) { return e1.weight < e2.weight; });
*/
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
    uint64_t bucketCnt; // = std::stoull(argv[3]);
    uint64_t numSamples;
    std::string eqClsListFile;
};

void loadEqs(std::string filename, std::vector<sdsl::rrr_vector < 63>> &bvs) {
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
    //BitVectorRRR bv(eqfile);
    std::cerr << "loaded all the equivalence classes: "
              << ((bvs.size() - 1) * EQS_PER_SLOT + bvs.back().size())
              << "\n";
}

void buildColor(std::vector<sdsl::rrr_vector < 63>>&bvs, std::vector<uint64_t> &eq,
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


int main(int argc, char *argv[]) {
    /* Let us create above shown weighted
       and undirected graph */
    using namespace clipp;
    enum class mode {
        build, ccInfo, help
    };
    mode selected = mode::help;
    Opts opt;
    auto ccInfo_mode = (
            command("ccInfo").set(selected, mode::ccInfo),
                    required("-e", "--edge-filename") &
                    value("edge_filename", opt.filename) % "file containing list of eq. class edges.",
                    required("-n", "--eqCls-cnt") &
                    value("equivalenceClass_count", opt.numNodes) % "Total number of equivalence (color) classes.",
                    required("-b", "--bucket-cnt") &
                    value("bucket_count", opt.bucketCnt) % "Total number of valid distances.",
                    required("-s", "--numSamples") &
                    value("numSamples", opt.numSamples) % "Total number of experiments (samples)."
    );

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
            (build_mode | ccInfo_mode | command("help").set(selected, mode::help)
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

    //explore_options_verbose(res);

    if (res) {
        switch (selected) {
            //case mode::ccInfo: query_main(qopt);  break;
            //case mode::build: validate_main(vopt);  break;
            case mode::help:
                std::cerr << make_man_page(cli, "MSF");
                break;
        }
    }

    std::cerr << "here are the inputs: \n"
              << opt.filename << "\n"
              << opt.numNodes << "\n"
              << opt.bucketCnt << "\n";
    opt.numNodes++; // number of nodes is one more than input including the zero node
    ifstream file(opt.filename);

    Graph g(opt.bucketCnt);

    uint32_t w_;
    uint32_t n1, n2, edgeCntr{0}, zero{(uint32_t) opt.numNodes - 1};
    std::string tmp;
    {
        //unordered_set<uint64_t> nodes;
        sdsl::bit_vector nodes(opt.numNodes, 0);
        std::getline(file, tmp);
        for (uint64_t i = 0; i < opt.numNodes; i++) {
            g.addEdge(zero, i, opt.numSamples);
        }
        while (file.good()) {
            file >> n1 >> n2 >> w_;
            g.addEdge(n1, n2, w_);
            nodes[n1] = 1;
            nodes[n2] = 1;
            edgeCntr++;
            //nodes.insert(n1);
            //nodes.insert(n2);
        }
        //file.clear();
        //file.seekg(0, ios::beg);
        file.close();

        uint64_t distinctNodes{0};
        for (uint64_t i = 0; i < nodes.size(); i += 64) {
            distinctNodes += sdsl::bits::cnt(nodes.get_int(i, 64));
        }
        std::cerr << "# of nodes: " << distinctNodes//nodes.size()
                  << "\n# of edges: " << edgeCntr
                  << "\n";
//        nodes.clear();
    }
    g.V = opt.numNodes;
    g.mst.resize(opt.numNodes);
    g.E = edgeCntr;
    //g.V = numNodes;
    //ifstream file(filename);

    DisjointSets ds = g.kruskalMSF(opt.bucketCnt/*, file*/);

    if (selected == mode::ccInfo) {
        for (auto &el : ds.els) {
            if (el.w != 0) {
                std::cout << el.edges << "\t" << el.w << "\t" << el.rnk << "\n";
            }
        }
    }

    if (selected == mode::build) {
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
                // just to confirm this, I'll continue
                if (check) {
                    std::cerr << "A TERRIBLE BUG!!!\n"
                              << "finding a node with no edges should only happen once at the root\n";
                    std::exit(1);
                }
                check = true;
                root = node;
                continue;
            }
            auto buck = g.mst[node][0].bucket;
            auto idx = g.mst[node][0].idx;
            uint64_t dest = g.edges[buck][idx].n1;
            uint64_t src = g.edges[buck][idx].n2;
            if (src == node) { // swap src & dest since src is the leaf
                std::swap(src, dest);
            }
            p2c[src].emplace_back(dest, g.edges[buck][idx].weight);
            g.mst[dest].erase(std::remove_if(
                    g.mst[dest].begin(), g.mst[dest].end(),
                    [buck, idx](const EdgePtr &x) {
                        return x.bucket == buck && x.idx == idx;
                    }), g.mst[dest].end());
            g.mst[src].erase(std::remove_if(
                    g.mst[src].begin(), g.mst[src].end(),
                    [buck, idx](const EdgePtr &x) {
                        return x.bucket == buck && x.idx == idx;
                    }), g.mst[src].end());
            if (g.mst[src].size() == 1) {
                q.push(src);
            }
            nodeCntr++;
            if (nodeCntr % 10000000 == 0) {
                std::cerr << nodeCntr << " nodes processed toward root\n";
            }
        }

        // calculate all the stats!!
        // create the data structures
        std::cerr << "Calculate Stats .. \n";
        nodeCntr = 0;
        std::queue<Path> pq;
        pq.push(Path(root, 0, 0));
        double avgDegree{0};
        uint64_t internalNodeCnt{0};
        std::cout << "-1\t" << root << "\t" << 0 << "\t" << 0 << "\n";
        while (!pq.empty()) {
            Path p = pq.front();
            pq.pop();
            avgDegree += p2c[p.id].size();
            if (p2c[p.id].size() != 0) {
                internalNodeCnt++;
            }
            for (auto &c : p2c[p.id]) {
                Path cp(c.id, p.steps + 1, p.weight + c.weight);
                std::cout << p.id << "\t" << cp.id << "\t" << cp.steps << "\t" << cp.weight << "\n";
                pq.push(cp);
            }
            nodeCntr++;
            if (nodeCntr % 10000000 == 0) {
                std::cerr << nodeCntr << " nodes processed from root\n";
            }
        }

        std::cerr << "internal node count: " << internalNodeCnt << "\n"
                  << "total out degree: " << avgDegree << "\n"
                  << "average degree: " << avgDegree/internalNodeCnt << "\t" << avgDegree/opt.numNodes << "\n";

    }

    return 0;
}

