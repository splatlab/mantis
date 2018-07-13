//
// Created by Fatemeh Almodaresi on 7/11/18.
//
// The algorithm's basic implementation taken from
// https://www.geeksforgeeks.org/kruskals-minimum-spanning-tree-using-stl-in-c/
//
#include<bits/stdc++.h>
#include <sstream>
#include <unordered_set>

using namespace std;

struct Edge {
    uint64_t n1;
    uint64_t n2;
    uint32_t weight;

    Edge(uint64_t inN1, uint64_t inN2, uint32_t inWeight)
            : n1(inN1), n2(inN2), weight(inWeight) {}
};

struct DisjointSetNode {
    uint64_t parent{0}, rnk{0}, w{0}, edges{0};

    void setParent(uint64_t p) { parent = p; }

    void mergeWith(DisjointSetNode &n, uint16_t edgeW) {
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
    int find(uint64_t u) {
        /* Make the parent of the nodes in the path
           from u--> parent[u] point to parent[u] */
        if (u != els[u].parent)
            els[u].parent = find(els[u].parent);
        return els[u].parent;
    }

    // Union by rank
    void merge(uint64_t x, uint64_t y, uint16_t edgeW) {
        x = find(x), y = find(y);

        /* Make tree with smaller height
           a subtree of the other tree  */
        if (els[x].rnk > els[y].rnk) {
            els[x].mergeWith(els[y], edgeW);

        } else {// If rnk[x] <= rnk[y]
            els[y].mergeWith(els[x], edgeW);
        }
    }
};

// Structure to represent a graph
struct Graph {
    uint64_t V, E;
    vector<Edge> edges;

    Graph() {}

    // Utility function to add an edge
    void addEdge(uint64_t u, uint64_t v, uint16_t w) {
        edges.emplace_back(u, v, w);
    }

    // Function to find MST using Kruskal's
    // MST algorithm
    DisjointSets kruskalMSF(uint64_t bucketCnt, ifstream &file) {
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
        // Iterate through all sorted edges
        for (auto bucketCntr = 1; bucketCntr <= bucketCnt; bucketCntr++) {
            //ifstream file(filename);
            std::getline(file, tmp);
            while (file.good()) {
                file >> n1 >> n2 >> w;
                if (w == bucketCntr) {
                    uint64_t u = n1;
                    uint64_t v = n2;
                    uint64_t set_u = ds.find(u);
                    uint64_t set_v = ds.find(v);

                    // Check if the selected edge is creating
                    // a cycle or not (Cycle is created if u
                    // and v belong to same set)
                    if (set_u != set_v) {
                        // Current edge will be in the MST
                        // Merge two sets
                        ds.merge(set_u, set_v, w);
                        mergeCntr++;
                    }/* else {
                            if (nodes.find(u) == nodes.end() || nodes.find(v) == nodes.end())
                                std::cerr << u << " " << v << " " << set_u << " " << set_v << "\n";
                    }*/
                    cntr++;
                    if (cntr % 1000000 == 0) {
                        std::cerr << "edge " << cntr << " " << mergeCntr << "\n";
                    }
                }
            }
            file.clear();
            file.seekg(0, file.beg);

        }
        file.close();
        std::cerr << "final # of edges: " << cntr
                  << " # of merges: " << mergeCntr << "\n";
        return ds;
    }
};

int main(int argc, char *argv[]) {
    /* Let us create above shown weighted
       and undirected graph */
    std::string filename = argv[1];
    uint64_t numNodes = std::stoull(argv[2]);
    uint64_t bucketCnt = std::stoull(argv[3]);

    ifstream file(filename);

    Graph g;

    /*uint32_t w_;
    uint64_t n1, n2;
    std::string tmp;
    unordered_set<uint64_t> nodes;
    std::getline(file, tmp);
    while (file.good()) {
        file >> n1 >> n2 >> w_;
        g.addEdge(n1, n2, w_);
        nodes.insert(n1);
        nodes.insert(n2);
    }
    //file.clear();
    //file.seekg(0, ios::beg);
    file.close();*/
    /* std::cerr << "# of nodes: " << nodes.size()
               << "\n# of edges: " << g.edges.size()
               << "\n";
     g.V = numNodes;
     g.E = g.edges.size();
     nodes.clear();*/
    g.V = numNodes;
    //ifstream file(filename);

    DisjointSets ds = g.kruskalMSF(bucketCnt, file);

    for (auto &el : ds.els) {
        if (el.w != 0) {
            std::cout << el.edges << "\t" << el.w << "\n";
        }
    }
    return 0;
}

