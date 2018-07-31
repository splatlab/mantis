//
// Created by Fatemeh Almodaresi on 7/11/18.
//
// The algorithm's basic implementation taken from
// https://www.geeksforgeeks.org/kruskals-minimum-spanning-tree-using-stl-in-c/
//

#include <set>
#include <queue>
#include "MSF.h"
#include "sparsepp/spp.h"

struct Hop {
    uint32_t id;
    uint32_t dist;
    uint16_t level;

    Hop(uint32_t id, uint32_t dist, uint16_t level) : id(id), dist(dist), level(level) {}

    bool operator<(const Hop &rhs) const {
        return level < rhs.level;
    }

};

struct Opts {
    std::string filename;
    uint64_t numNodes; // = std::stoull(argv[2]);
    uint16_t numSamples;
    std::string eqClsListFile;
    std::string outputDir;
    uint16_t hops;
};


int main(int argc, char *argv[]) {
    /* Let us create above shown weighted
       and undirected graph */
    using namespace clipp;
    enum class mode {
        build, fillGraph, help
    };
    mode selected = mode::help;
    Opts opt;

    auto build_mode = (
            command("build").set(selected, mode::build),
                    required("-e", "--edge-filename") &
                    value("edge_filename", opt.filename) % "File containing list of eq. class edges.",
                    required("-n", "--eqCls-cnt") &
                    value("equivalenceClass_count", opt.numNodes) % "Total number of equivalence (color) classes.",
                    required("-s", "--numSamples") &
                    value("numSamples", opt.numSamples) % "Total number of experiments (samples).",
                    required("-c", "--eqCls-lst") &
                    value("eqCls_list_filename", opt.eqClsListFile) %
                    "File containing list of equivalence (color) classes.",
                    required("-o", "--output_dir") &
                    value("output directory", opt.outputDir) % "Directory that all the int_vectors will be stored in."
    );

    auto fillGraph_mode = (
            command("fillGraph").set(selected, mode::fillGraph),
                    required("-e", "--edge-filename") &
                    value("edge_filename", opt.filename) % "File containing list of eq. class edges.",
                    required("-s", "--numSamples") &
                    value("numSamples", opt.numSamples) % "Total number of experiments (samples).",
                    required("-n", "--eqCls-cnt") &
                    value("equivalenceClass_count", opt.numNodes) % "Total number of equivalence (color) classes.",
                    required("-c", "--eqCls-lst") &
                    value("eqCls_list_filename", opt.eqClsListFile) %
                    "File containing list of equivalence (color) classes.",
                    required("-o", "--output_dir") &
                    value("output directory", opt.outputDir) % "Directory that all the int_vectors will be stored in.",
                    required("-h", "--hops") &
                    value("hops", opt.hops) % "# of hops to search for each node for a direct link with smaller weight."
    );

    auto cli = (
            (build_mode | fillGraph_mode | command("help").set(selected, mode::help)
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
              << opt.numSamples << "\n";
    uint64_t numWrds = (uint64_t) std::ceil((double) opt.numSamples / 64.0);
    opt.numNodes++; // number of nodes is one more than input including the zero node

    std::cerr << "Loading all the equivalence classes first .. \n";
    std::vector<sdsl::rrr_vector < 63>>
    eqs;
    eqs.reserve(20);
    loadEqs(opt.eqClsListFile, eqs);
    std::cerr << "Done loading list of equivalence class buckets\n";

    ifstream file(opt.filename);

    uint16_t w_;
    uint32_t n1, n2, edgeCntr{0}, zero{(uint32_t) opt.numNodes - 1};
    if (selected == mode::build) {
        Graph g(opt.numSamples);

        {
            std::cerr << "Adding edges from 0 to each node with number of set bits in the node as weight .. \n";
            for (uint32_t i = 0; i < opt.numNodes; i++) {
                uint16_t ones = sum1s(eqs, i, opt.numSamples, numWrds);
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

        DisjointSets ds = g.kruskalMSF(opt.numSamples);
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
        std::cerr << "Creating parentBV...\n";
        sdsl::int_vector<> parentbv(opt.numNodes, 0, ceil(log2(opt.numNodes)));

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
                parentbv[node] = node;
                check = true;
                // Update the total weight to contain the delta of root from bv of zero
                uint16_t ones = 1; // If the root is zero itself
                if (node != zero) // otherwise
                    ones = sum1s(eqs, node, opt.numSamples, numWrds);
                g.mst_totalWeight += ones;
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
            if (dest == zero) { // we break the tree from node zero
                parentbv[zero] = zero;
                // we're breaking the tree, so should remove the edge weight from total weights
                // But instead we're gonna spend one slot to store 0 as the delta of the zero node from zero node
                g.mst_totalWeight = g.mst_totalWeight - g.edges[buck][idx].weight + 1;
            } else {
                parentbv[dest] = src;
            }
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
            // but if the src has turned to a leaf (node with one edge)
            // add it to the queue
            if (g.mst[src].size() == 1) {
                q.push(src);
            }
            nodeCntr++; // just a counter for the log
            if (nodeCntr % 10000000 == 0) {
                std::cerr << nodeCntr << " nodes processed toward root\n";
            }
        }

        // create the data structures
        std::cerr << "Sum of MST weights: " << g.mst_totalWeight << "\n";
        std::cerr << "Creating deltaBV and bBV...\n";
        sdsl::int_vector<> deltabv(g.mst_totalWeight, 0, ceil(log2(opt.numSamples)));
        sdsl::bit_vector bbv(g.mst_totalWeight, 0);

        uint64_t deltaOffset{0};
        for (uint64_t i = 0; i < parentbv.size(); i++) {
            std::vector<uint32_t> deltas;
            if (i == zero) {
                deltaOffset++;
            } else if (parentbv[i] == zero || parentbv[i] == i) {
                deltas = getDeltaList(eqs, i, opt.numSamples, numWrds);
            } else {
                deltas = getDeltaList(eqs, parentbv[i], i, opt.numSamples, numWrds);
            }
            for (auto &v : deltas) {
                deltabv[deltaOffset] = v;
                deltaOffset++;
            }
            bbv[deltaOffset - 1] = 1;
            if (i % 10000000 == 0) {
                std::cerr << i << " nodes parents processed\n";
            }
        }


        std::cerr << "Sum of MST weights: " << g.mst_totalWeight << "\n";

        sdsl::store_to_file(parentbv, opt.outputDir + "/parents.bv");
        sdsl::store_to_file(deltabv, opt.outputDir + "/deltas.bv");
        sdsl::store_to_file(bbv, opt.outputDir + "/boundary.bv");
    } else if (selected == mode::fillGraph) {
        //spp::sparse_hash_map<uint32_t, std::vector<std::pair<uint32_t, uint16_t>>> nodes;
        std::vector<std::vector<std::pair<uint32_t, uint16_t>>> nodes(opt.numNodes);
        std::cerr << "Adding edges between color classes .. \n";
        while (file.good()) {
            file >> n1 >> n2 >> w_;
            nodes[n1].push_back(std::make_pair(n2, w_));
            nodes[n2].push_back(std::make_pair(n1, w_));
            edgeCntr++;
            if (edgeCntr % 10000000 == 0) {
                std::cerr << edgeCntr << " passed\n";
            }
        }
        file.close();
        std::ofstream of(opt.outputDir + "/extraEdges.lst");
        std::cerr << "Done adding edges between color classes .. \n";
        std::cerr << "Total # of nodes : " << nodes.size() << "\n";
        std::cerr << "Hops required: " << opt.hops << "\n";
        uint32_t nodeCntr{0};
        uint64_t eqWrds{static_cast<uint64_t>(((opt.numSamples - 1) / 64) + 1)};
        //std::vector<vector<uint32_t>> nodeEqs(nodes.size());
        std::cerr << "initialized\n";
        for (auto &neis : nodes) {
            std::unordered_set<uint32_t> visited;
            visited.insert(nodeCntr);
            std::priority_queue<Hop> hops;
            std::vector<uint64_t> eq1(eqWrds);
            buildColor(eqs, eq1, nodeCntr, opt.numSamples);

            //std::cerr << "n " << nodeCntr << " " << neis.size() << "\n";
            uint64_t maxNei{100};
            if (neis.size() > maxNei) continue;
            for (auto &nei : neis) {
                visited.insert(nei.first);
                Hop nh(nei.first, nei.second, 1);
                hops.push(nh);
            }

            while (!hops.empty() && hops.top().level != opt.hops) {
                auto nei = hops.top();
                hops.pop();
                //std::cerr << nei.level << " " << nei.id << " " << nei.dist << " " << nodes[nei.id].size() << "\n";
                if (nodes[nei.id].size() > maxNei) continue;
                for (auto &neinei : nodes[nei.id]) {
                    //std::cerr << "  " << neinei.first << "\n";
                    if (visited.find(neinei.first) == visited.end()) {
                        visited.insert(neinei.first);
                        Hop nh(neinei.first, neinei.second + nei.dist, nei.level + 1);
                        if (nodeCntr < nh.id) {
                            std::vector<uint64_t> eq2(eqWrds);
                            buildColor(eqs, eq2, nh.id, opt.numSamples);
                            uint64_t directDist = 0;//hammingDist(eqs, id, nh.id, opt.numSamples);
                            for (uint64_t i = 0; i < eq1.size(); i++) {
                                if (eq1[i] != eq2[i])
                                    directDist += sdsl::bits::cnt(eq1[i] ^ eq2[i]);
                            }
                            if (directDist < nh.dist) {
                                of << nodeCntr << "\t" << nh.id << "\t" << directDist << "\n";
                            }
                        }
                        hops.push(nh);
                    }
                }
            }
            nodeCntr++;
            if (nodeCntr % 1000000 == 0) {
                std::cerr << nodeCntr << " passed\n";
            }
        }
    }
    return 0;
}

