#include <utility>

//
// Created by Fatemeh Almodaresi on 2018-10-13.
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

void MST::buildMST() {
    buildEdgeSets();
    calculateWeights();
    findMST();
}

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

void MST::findNeighborEdges(CQF<KeyObject>& cqf, CQF<KeyObject>::Iterator it) {
    KeyObject keyobj = *it;
    //node curn(k, hash_64i(keyobj.key, BITMASK(cqf->keybits())));
    dna::canonical_kmer curr_node(static_cast<int>(k), keyobj.key);
    workItem cur = {curr_node, static_cast<uint32_t >(keyobj.count - 1)};
    uint64_t neighborCnt{0};
    for (auto &nei : neighbors(cqf, cur)) {
        neighborCnt++;
        if (cur.colorId < nei.colorId) {
            Edge e(static_cast<uint32_t>(cur.colorId), static_cast<uint32_t>(nei.colorId));

            auto& edgeset = edgesetList[fetchBucketId(cur.colorId, nei.colorId)];
            if (edgeset.find(e) == edgeset.end()) {
                edgeset.insert(e);
                num_edges++;
            }
        }
    }
}

bool MST::calculateWeights() {
    std::vector<std::string> eqclass_files =
            mantis::fs::GetFilesExt(prefix.c_str(), mantis::EQCLASS_FILE);
    return true;
}

bool MST::findMST() {
    return true;
}

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

bool MST::exists(CQF<KeyObject>& cqf, dna::canonical_kmer e, uint64_t &eqid) {
    KeyObject key(e.val, 0, 0);
    auto eqidtmp = cqf.query(key, 0 /*QF_KEY_IS_HASH | QF_NO_LOCK*/);
    if (eqid) {
        eqid = eqidtmp - 1;
        return true;
    }
    return false;
}

// c1 <= c2
inline uint64_t MST::fetchBucketId(uint64_t c1, uint64_t c2) {
    uint64_t cb1 = c1/mantis::NUM_BV_BUFFER;
    uint64_t cb2 = c2/mantis::NUM_BV_BUFFER;
    return cb1*num_of_ccBuffers+cb2;
}

int build_mst_main( QueryOpts& opt ) {
    MST mst(opt.prefix);
    mst.buildMST();
    return 0;
}