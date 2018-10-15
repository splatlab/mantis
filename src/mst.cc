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

    eqclass_files =
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
    encodeColorClassUsingMST();
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
    sdsl::bit_vector nodes(num_of_ccBuffers*mantis::NUM_BV_BUFFER, 0);
    auto it = cqf.begin();
    while (!it.done()) {
        KeyObject keyObject = *it;
        nodes[keyObject.count-1] = 1; // set the seen color class id bit
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

    // Add an edge between edch color class ID and node zero
    console->info("Adding edges from dummy node zero to each color class Id");
    zero = static_cast<colorIdType>(num_colorClasses);
    num_colorClasses++; // zero is now a dummy color class with ID equal to actual num of color classes
    for (colorIdType colorId = 0; colorId < num_colorClasses; colorId++) {
        edgesetList[getBucketId(colorId, zero)].insert(Edge(colorId, zero));
    }
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

    console->info("Going over all the edges and calculating the weights.");
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
                auto w = hammingDist(edge.n1, edge.n2); // hammingDist uses bv1 and bv2
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
            colorIdType u = it->n1;
            colorIdType v = it->n2;
            colorIdType root_of_u = ds.find(u);
            colorIdType root_of_v = ds.find(v);

            // Check if the selected edge is causing a cycle or not
            // (A cycle is induced if u and v belong to the same set)
            if (root_of_u != root_of_v) {
                // Merge two sets
                ds.merge(root_of_u, root_of_v, w);
                // Current edge will be in the MST
                mst[u].emplace_back(v, w);
                mst[v].emplace_back(u, w);
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
    mstTotalWeight++;//1 empty slot for root (zero)
    console->info("MST Construction finished:"
                  "\n\t# of edges: {}"
                  "\n\t# of merges: {}"
                  "\n\tmst weight sum: {}",
                  edgeCntr, selectedEdgeCntr, mstTotalWeight);
    return ds;
}

bool MST::encodeColorClassUsingMST() {
    // build mst of color class graph
    kruskalMSF();

    uint64_t nodeCntr{0};
    // encode the color classes using mst
    std::cerr << "Creating parentBV...\n";
    sdsl::int_vector<> parentbv(num_colorClasses, 0, ceil(log2(num_colorClasses)));
    // create and fill the deltabv and boundarybv data structures
    sdsl::bit_vector bbv;
    {// putting weightbv inside the scope so its memory is freed after we're done with it
        sdsl::int_vector<> weightbv(num_colorClasses, 0, ceil(log2(numSamples)));

        {
            sdsl::bit_vector visited(num_colorClasses, 0);
            bool check = false;
            std::queue<colorIdType> q;
            q.push(zero); // Root of the tree is zero
            parentbv[zero] = zero; // and it's its own parent (has no parent)
            while (!q.empty()) {
                colorIdType parent = q.front();
                q.pop();
                for (auto &neighbor :mst[parent]) {
                    if (!visited[neighbor.first]) {
                        parentbv[neighbor.first] = parent;
                        weightbv[neighbor.first] = neighbor.second;
                        q.push(neighbor.first);
                    }
                }
                visited[parent] = 1;
                nodeCntr++; // just a counter for the log
                if (nodeCntr % 10000000 == 0) {
                    std::cerr << "\rset parent of " << nodeCntr << " ccs";
                }
            }
        }

        // filling bbv
        // resize bbv
        nodeCntr=0;
        bbv.resize(mstTotalWeight);
        uint64_t deltaOffset{0};
        for (uint64_t i = 0; i < num_colorClasses; i++) {
            std::vector<uint32_t> deltas;
            if (i == parentbv[i]) { //it's the root (zero here)
                deltaOffset++;
            } else {
                deltaOffset += weightbv[i];
            }
            bbv[deltaOffset - 1] = 1;
            if (i % 10000000 == 0) {
                std::cerr << "\rset delta vals for " << nodeCntr << " ccs";
            }
        }
    }
    sdsl::int_vector<> deltabv(mstTotalWeight, 0, ceil(log2(numSamples)));
    sdsl::bit_vector visited(num_colorClasses, 0);

    // fill in deltabv
    sdsl::bit_vector::select_1_type sbbv = sdsl::bit_vector::select_1_type(&bbv);
    console->info("Start going over all the edges and calculating the weights.");
    for (auto i = 0; i < eqclass_files.size(); i++) {
        sdsl::load_from_file(*bv1, eqclass_files[i]);
        for (auto j = i; j < eqclass_files.size(); j++) {
            std::cerr << "\rset delta vals for cc buffers " << i << " & " << j;
            if (i == j) {
                bv2 = bv1;
            } else {
                sdsl::load_from_file(*bv2, eqclass_files[j]);
            }
            for (colorIdType p=0; p < parentbv.size(); p++) {
                if (getBucketId(p, parentbv[p]) == i*mantis::NUM_BV_BUFFER+j) {
                    auto deltaOffset = (p > 0) ? (sbbv(p) + 1) : 0;
                    for (auto &v : getDeltaList(p, parentbv[p])) {
                        deltabv[deltaOffset] = v;
                        deltaOffset++;
                    }
                }
            }
        }
    }
    delete bv1;
    delete bv2;

    sdsl::store_to_file(parentbv, std::string(prefix + mantis::PARENTBV_FILE));
    sdsl::store_to_file(deltabv, std::string(prefix + mantis::DELTABV_FILE) );
    sdsl::store_to_file(bbv, std::string(prefix + mantis::BOUNDARYBV_FILE) );
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
    workItem cur = {curr_node, static_cast<colorIdType>(keyobj.count - 1)};
    uint64_t neighborCnt{0};
    for (auto &nei : neighbors(cqf, cur)) {
        neighborCnt++;
        if (cur.colorId < nei.colorId) {
            Edge e(static_cast<colorIdType>(cur.colorId), static_cast<colorIdType>(nei.colorId));

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

//
/**
 * for two non-zero nodes, list indices that xor of the bits is 1
 * for one non-zero node, list indices that the bit is 1
 *
 * @param eqid1
 * @param eqid2
 * @return delta list
 */
std::vector<uint32_t> MST::getDeltaList(uint64_t eqid1,uint64_t eqid2) {
    std::vector<uint32_t> res;
    std::vector<uint64_t> eq1(((numSamples - 1) / 64) + 1, 0), eq2(((numSamples - 1) / 64) + 1, 0);
    buildColor(eq1, eqid1, bv1);
    buildColor(eq2, eqid2, bv2);

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

/**
 * Loads the bitvector corresponding to eqId
 * @param eq list of words each representing 64 bits of eqId bv (output)
 * @param eqid color id
 * @param bv the large bv collapsing all eq ids color bv in a bucket
 */
void MST::buildColor(std::vector<uint64_t> &eq, uint64_t eqid, BitVectorRRR *bv) {
    if (eqid == zero) return;
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
    if (c1 > c2)
        std::swap(c1, c2);
    uint64_t cb1 = c1 / mantis::NUM_BV_BUFFER;
    uint64_t cb2 = c2 / mantis::NUM_BV_BUFFER;
    return cb1 * num_of_ccBuffers + cb2;
}

/**
 ********* MAIN *********
 * main function to call Color graph and MST construction and color class encoding and serializing
 */
int build_mst_main(QueryOpts &opt) {
    MST mst(opt.prefix);
    mst.buildMST();
    return 0;
}