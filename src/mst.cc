#include <utility>

//
// Created by Fatemeh Almodaresi.
//

#include <MantisFS.h>
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
    uint64_t kmerCntr{0};
    auto it = cqf.begin();
    while (!it.done()) {
        findNeighborEdges(cqf, it);
        ++it;
        kmerCntr++;
        if (kmerCntr % 10000000 == 0) {
            std::cerr << "\r" << kmerCntr << " kmers & " << num_edges << " edges\n";
        }
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
    std::vector<std::string> eqclass_files =
            mantis::fs::GetFilesExt(prefix.c_str(), mantis::EQCLASS_FILE);
    return true;
}

bool MST::findMST() {
    return true;
}

/**
 * for each element of cqf, finds its neighbors
 * and adds an edge of the element's colorId and its neighbor's
 * @param cqf (required to query for existence of neighbors)
 * @param it iterator to the elements of cqf
 */
void MST::findNeighborEdges(CQF<KeyObject>& cqf, CQF<KeyObject>::Iterator it) {
    KeyObject keyobj = *it;
    dna::canonical_kmer curr_node(static_cast<int>(k), keyobj.key);
    workItem cur = {curr_node, static_cast<uint32_t >(keyobj.count - 1)};
    uint64_t neighborCnt{0};
    for (auto &nei : neighbors(cqf, cur)) {
        neighborCnt++;
        if (cur.colorId < nei.colorId) {
            Edge e(static_cast<uint32_t>(cur.colorId), static_cast<uint32_t>(nei.colorId));

            auto& edgeset = edgesetList[getBucketId(cur.colorId, nei.colorId)];
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
std::set<workItem> MST::neighbors(CQF<KeyObject>& cqf, workItem n) {
    std::set<workItem> result;
    for (const auto b : dna::bases) {
        uint64_t eqid = 0;
        if (exists(cqf, n.node << b,eqid))
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
bool MST::exists(CQF<KeyObject>& cqf, dna::canonical_kmer e, uint64_t &eqid) {
    KeyObject key(e.val, 0, 0);
    auto eqidtmp = cqf.query(key, 0 /*QF_KEY_IS_HASH | QF_NO_LOCK*/);
    if (eqid) {
        eqid = eqidtmp - 1;
        return true;
    }
    return false;
}

/**
 * calculates the edge corresponding bucket id c1 <= c2
 * @param c1 first colorId
 * @param c2 second colorId
 * @return bucket id
 */
inline uint64_t MST::getBucketId(uint64_t c1, uint64_t c2) {
    uint64_t cb1 = c1/mantis::NUM_BV_BUFFER;
    uint64_t cb2 = c2/mantis::NUM_BV_BUFFER;
    return cb1*num_of_ccBuffers+cb2;
}

int build_mst_main( QueryOpts& opt ) {
    MST mst(opt.prefix);
    mst.buildMST();
    return 0;
}