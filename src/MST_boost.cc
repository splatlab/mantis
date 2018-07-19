//
// Created by Fatemeh Almodaresi on 7/19/18.
//


#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/kruskal_min_spanning_tree.hpp>
#include <boost/range/irange.hpp>
#include <iostream>
#include <fstream>
#include "clipp.h"

struct Opts {
    std::string filename;
    uint64_t numNodes; // = std::stoull(argv[2]);
    uint64_t bucketCnt; // = std::stoull(argv[3]);
    uint64_t numSamples;
    std::string eqClsListFile;
};

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
                    value("bucket_count", opt.bucketCnt) % "Total number of valid distances."
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

    using namespace boost;
    typedef adjacency_list <vecS, vecS, undirectedS, no_property, property<edge_weight_t, uint32_t>> Graph;
    typedef graph_traits<Graph>::edge_descriptor Edge;
    typedef graph_traits<Graph>::vertex_descriptor Vertex;
    typedef std::pair<uint64_t, uint64_t> E;

    if (selected == mode::build) {

        Graph g(opt.numNodes+1);
        property_map<Graph, edge_weight_t>::type weightmap = get(edge_weight, g);
        std::ifstream file(opt.filename);
        uint64_t n1,n2, zero{opt.numNodes+2}, edgeCntr{0};
        uint32_t w;
        for (auto i : irange((uint64_t)0, opt.numNodes)) {
            add_edge(0, i, opt.numSamples, g);
            edgeCntr++;
            if (edgeCntr % 1000000 == 0) {

            }
        }
        std::cerr << edgeCntr << " zero-end edges\n";
        edgeCntr = 0;
        while (file.good()) {
            file >> n1 >> n2 >> w;
            /*Edge e;
            bool inserted;*/
            /*tie(e, inserted) = */add_edge(n1, n2, w, g);
            //weightmap[e] = w;
            edgeCntr++;
            if (edgeCntr % 10000000 == 0) {
                std::cerr << edgeCntr << "\n";
            }
        }
        std::cerr << edgeCntr << " edges\n";

        std::vector<Edge> spanning_tree;

        kruskal_minimum_spanning_tree(g, std::back_inserter(spanning_tree));

        std::cerr << "Done building MST\nTotal # of edges: " << spanning_tree.size() << "\n";
        /*std::cout << "Print the edges in the MST:" << std::endl;
        for (std::vector<Edge>::iterator ei = spanning_tree.begin();
             ei != spanning_tree.end(); ++ei) {
            std::cout << source(*ei, g) << " <--> " << target(*ei, g)
                      << " with weight of " << weight[*ei]
                      << std::endl;
        }*/
    }
    return EXIT_SUCCESS;
}