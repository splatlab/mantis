#include <utility>

//
// Created by Fatemeh Almodaresi.
//
#include <string>
#include <sstream>

#include "MantisFS.h"
#include "mst.h"
#include "ProgOpts.h"

MST::MST(std::string prefixIn) : prefix(std::move(prefixIn)) {

    console = spdlog::stdout_color_mt("mantis_console");

    // Make sure the prefix is a full folder
    if (prefix.back() != '/') {
        prefix.push_back('/');
    }
    // make the output directory if it doesn't exist
    if (!mantis::fs::DirExists(prefix.c_str())) {
        console->error("Index parent directory {} does not exist", prefix);
        std::exit(1);
    }

    std::vector<std::string> eqclass_files =
            mantis::fs::GetFilesExt(prefix.c_str(), mantis::EQCLASS_FILE);
    num_of_ccBuffers = eqclass_files.size();

    std::string sample_file(prefix.c_str(), mantis::SAMPLEID_FILE);
    std::ifstream sampleid(sample_file.c_str());
    uint32_t tmp;
    while (sampleid >> tmp >> tmp) {
        numSamples++;
    }
    sampleid.close();
}

/**
 * Building an MST consists of 3 main steps:
 * 1. construct the color graph for all the colorIds derived from dbg
 *      This phase just requires loading the CQF
 * 2. calculate the weights of edges in the color graph
 *      This phase requires at most two buffers of color classes
 * 3. find MST of the weighted color graph
 */
void MST::buildMST() {
    buildEdgeSets();
    calculateWeights();
    findMST();
}

/**
 * iterate over all elements of CQF,
 * find all the existing neighbors, and build a color graph based on that
 * @return true if the color graph build was successful
 */
bool MST::buildEdgeSets() {
    console->info("Reading colored dbg from disk.");
    std::string cqf_file(prefix + mantis::CQF_FILE);
    CQF<KeyObject> cqf(cqf_file, CQF_FREAD);
    k = cqf.keybits() / 2;
    console->info("Done Reading colored dbg from disk. k is {}", k);
    console->info("Starting to iterate over cqf ...");
    uint64_t kmerCntr{0};
    zero = num_of_ccBuffers*mantis::NUM_BV_BUFFER; // assuming last buffer is full
    sdsl::bit_vector nodes(zero, 0);
    auto it = cqf.begin();
    while (!it.done()) {
        KeyObject keyObject = *it;
        // Add an edge between the color class and node zero
        uint32_t c1{static_cast<uint32_t >(keyObject.count-1)}, c2{static_cast<uint32_t >(zero)};
        if (nodes[c1] == 0) { // haven't seen the color class ID yet
            nodes[c1] = 1;
            if (c1 == c2) {
                console->error("We shouldn't have a color class with ID {}", zero);
                std::exit(1);
            }
            if (c1 > c2) {
                std::swap(c1, c2);
            }
            edgesetList[getBucketId(c1, c2)].insert(Edge(c1, c2));
        }
        // Add an edge between the color class and each of its neighbors' colors in dbg
        findNeighborEdges(cqf, keyObject);
        ++it;
        kmerCntr++;
        if (kmerCntr % 10000000 == 0) {
            std::cerr << "\r" << kmerCntr << " kmers & " << num_edges << " edges\n";
        }
    }

    // count total number of color classes:
    for (uint64_t i = 0; i < nodes.size(); i += 64) {
        num_colorClasses += sdsl::bits::cnt(nodes.get_int(i, 64));
    }
    zero = num_colorClasses;
    return true;
}

/**
 * load the color class table in parts
 * calculate the hamming distance between the color bitvectors fetched from color class table
 * for each pair of color IDs
 * having w buckets where w is the maximum possible weight (number of experiments)
 * put the pair in its corresponding bucket based on the hamming distance value (weight)
 * @return true if successful
 */
bool MST::calculateWeights() {
    std::vector<std::string> eqclass_files =
            mantis::fs::GetFilesExt(prefix.c_str(), mantis::EQCLASS_FILE);

    // sort eqclass_files
    // note to @robP: It terribly statically relies on the format of the input files!!
    std::sort(eqclass_files.begin(), eqclass_files.end(), [this](std::string &s1, std::string &s2) {
        uint32_t id1, id2;
        std::stringstream ss1(first_part(last_part(s1, '/'), '_'));
        std::stringstream ss2(first_part(last_part(s2, '/'), '_'));
        if ((ss1 >> id1).fail() || !(ss1 >> std::ws).eof() ||
            (ss2 >> id2).fail() || !(ss2 >> std::ws).eof()) {
            console->error("file name does not start with a number : {}, {}", s1, s2);
        }
        return id1 < id2;
    });

    console->info("Start going over all the edges and calculating the weights.");
    for (auto i = 0; i < eqclass_files.size(); i++) {
        sdsl::load_from_file(*bv1, eqclass_files[i]);
        for (auto j = i; j < eqclass_files.size(); j++) {
            auto &edgeset = edgesetList[i * num_of_ccBuffers + j];
            if (i == j) {
                bv2 = bv1;
            } else {
                sdsl::load_from_file(*bv2, eqclass_files[j]);
            }
            for (auto &edge : edgeset) {
                auto w = hammingDist(edge.n1, edge.n2);
                if (w == 0) {
                    console->error("Hamming distance of 0 between edges {} & {}", edge.n1, edge.n2);
                    std::exit(1);
                }
                weightBuckets[w - 1].push_back(edge);
            }
            edgeset.clear();
        }
    }
    delete bv1;
    delete bv2;
    edgesetList.clear();
    return true;
}

bool MST::findMST() {
    // build mst
    kruskalMSF();
    // encode the color classes using mst

    return true;
}

/**
 * for each element of cqf, finds its neighbors
 * and adds an edge of the element's colorId and its neighbor's
 * @param cqf (required to query for existence of neighbors)
 * @param it iterator to the elements of cqf
 */
void MST::findNeighborEdges(CQF<KeyObject> &cqf, KeyObject &keyobj) {
    dna::canonical_kmer curr_node(static_cast<int>(k), keyobj.key);
    workItem cur = {curr_node, static_cast<uint32_t >(keyobj.count - 1)};
    uint64_t neighborCnt{0};
    for (auto &nei : neighbors(cqf, cur)) {
        neighborCnt++;
        if (cur.colorId < nei.colorId) {
            Edge e(static_cast<uint32_t>(cur.colorId), static_cast<uint32_t>(nei.colorId));

            auto &edgeset = edgesetList[getBucketId(cur.colorId, nei.colorId)];
            if (edgeset.find(e) == edgeset.end()) {
                edgeset.insert(e);
                num_edges++;
            }
        }
    }
}

/**
 * Find neighbors of a node in cqf
 * @param cqf
 * @param n : work_item containing node and colorId (colorId will be filled)
 * @return set of neighbors for current node n and their colorIds
 */
std::set<workItem> MST::neighbors(CQF<KeyObject> &cqf, workItem n) {
    std::set<workItem> result;
    for (const auto b : dna::bases) {
        uint64_t eqid = 0;
        if (exists(cqf, n.node << b, eqid))
            if (eqid != n.colorId)
                result.insert(workItem(n.node << b, eqid));
    }
    return result;
}

/**
 * searches for a kmer in cqf and returns the correct colorId if found
 * which is cqf count value - 1
 * @param cqf
 * @param e : search canonical kmer
 * @param eqid : reference to eqid that'll be set
 * @return true if eqid is found
 */
bool MST::exists(CQF<KeyObject> &cqf, dna::canonical_kmer e, uint64_t &eqid) {
    KeyObject key(e.val, 0, 0);
    auto eqidtmp = cqf.query(key, 0 /*QF_KEY_IS_HASH | QF_NO_LOCK*/);
    if (eqid) {
        eqid = eqidtmp - 1;
        return true;
    }
    return false;
}

/**
 * calculate hamming distance between the bvs of two color class ids
 * @param eqid1 first color class id
 * @param eqid2 second color class id
 * @return
 */
uint64_t MST::hammingDist(uint64_t eqid1, uint64_t eqid2) {
    uint64_t dist{0};
    std::vector<uint64_t> eq1(((numSamples - 1) / 64) + 1, 0), eq2(((numSamples - 1) / 64) + 1, 0);
    buildColor(eq1, eqid1, bv1);
    buildColor(eq2, eqid2, bv2);

    for (uint64_t i = 0; i < eq1.size(); i++) {
        if (eq1[i] != eq2[i])
            dist += sdsl::bits::cnt(eq1[i] ^ eq2[i]);
    }
    return dist;
}

/**
 * Loads the bitvector corresponding to eqId
 * @param eq list of words each representing 64 bits of eqId bv (output)
 * @param eqid color id
 * @param bv the large bv collapsing all eq ids color bv in a bucket
 */
void MST::buildColor(std::vector<uint64_t> &eq, uint64_t eqid, BitVectorRRR *bv) {
    if (eqid == 0) return;
    uint64_t i{0}, bitcnt{0}, wrdcnt{0};
    uint64_t offset = eqid % mantis::NUM_BV_BUFFER;
    while (i < numSamples) {
        bitcnt = std::min(numSamples - i, (uint64_t) 64);
        uint64_t wrd = bv->get_int(offset * numSamples + i, bitcnt);
        eq[wrdcnt++] = wrd;
        i += bitcnt;
    }
}

/**
 * calculates the edge corresponding bucket id c1 <= c2
 * @param c1 first colorId
 * @param c2 second colorId
 * @return bucket id
 */
inline uint64_t MST::getBucketId(uint64_t c1, uint64_t c2) {
    uint64_t cb1 = c1 / mantis::NUM_BV_BUFFER;
    uint64_t cb2 = c2 / mantis::NUM_BV_BUFFER;
    return cb1 * num_of_ccBuffers + cb2;
}

/**
 * Finding Minimim Spanning Forest of color graph using Kruskal Algorithm
 *
 * The algorithm's basic implementation taken from
 * https://www.geeksforgeeks.org/kruskals-minimum-spanning-tree-using-stl-in-c/
 * @return List of connected components in the Minimum Spanning Forest
 */
DisjointSets MST::kruskalMSF() {
    uint64_t numNodes = numSamples + 1;//for the dummy node with all bits=0
    uint32_t bucketCnt = numSamples;
    // Create disjoint sets
    DisjointSets ds(numNodes);

    uint64_t edgeCntr{0}, selectedEdgeCntr{0};
    uint32_t w{0};
    // Iterate through all sorted edges
    for (uint32_t bucketCntr = 0; bucketCntr < bucketCnt; bucketCntr++) {
        uint32_t edgeIdxInBucket = 0;
        for (auto it = weightBuckets[bucketCntr].begin(); it != weightBuckets[bucketCntr].end(); it++) {
            w = bucketCntr + 1;
            uint32_t u = it->n1;
            uint32_t v = it->n2;
            uint32_t root_of_u = ds.find(u);
            uint32_t root_of_v = ds.find(v);

            // Check if the selected edge is creating
            // a cycle or not (Cycle is induced if u
            // and v belong to the same set)
            if (root_of_u != root_of_v) {
                // Merge two sets
                ds.merge(root_of_u, root_of_v, w);
                // Current edge will be in the MST
                mst[u].emplace_back(v);
                mst[v].emplace_back(u);
                mstTotalWeight += w;
                selectedEdgeCntr++;
            }
            edgeCntr++;
            if (edgeCntr % 1000000 == 0) {
                console->info("\r{} edges processed and {} were selected", edgeCntr, selectedEdgeCntr);
            }
            edgeIdxInBucket++;
        }
        weightBuckets[bucketCntr].clear();
    }

    console->info("MST Construction finished:"
                  "\n\t# of edges: {}"
                  "\n\t# of merges: {}",
                  edgeCntr, selectedEdgeCntr);
    return ds;
}


int build_mst_main(QueryOpts &opt) {
    MST mst(opt.prefix);
    mst.buildMST();
    return 0;
}