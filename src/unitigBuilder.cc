//
// Created by Fatemeh Almodaresi on 2019-05-10.
//

#include "unitigBuilder.h"

UnitigBuilder::UnitigBuilder(std::string prefixIn, std::shared_ptr<spdlog::logger> loggerIn, uint32_t numThreads) :
        prefix(std::move(prefixIn)), nThreads(numThreads) {
    logger = loggerIn.get();
    logger->info("Reading cqf from disk.");
    std::string cqf_file(prefix + mantis::CQF_FILE);
    CQF<KeyObject> cqf(cqf_file, CQF_FREAD);
    k = cqf.keybits() / 2;
    logger->info("Done loading cdbg. k is {}", k);
    logger->info("Iterating over cqf & building unitigs ...");
    buildUnitigs(cqf);
}

/**
 * returns a set of prefix bases for the existing preceding kmers
 * @param cqf
 * @param node the kmer to search its previous neighbors (not canonicalized)
 * @return set of bases to prepend to "node" to build the existing neighbors
 */
std::set<dna::base> UnitigBuilder::prevNeighbors(CQF<KeyObject> &cqf, dna::kmer& node) {
    std::set<dna::base> result;
    for (const auto b : dna::bases) {
        colorIdType eqid = 0;
        if (exists(cqf, b >> node, eqid)) {
            result.insert(b);
        }
    }
    return result;
}

/**
* returns a set of prefix bases for the existing following kmers
* @param cqf
* @param node the kmer to search its succeeding neighbors (not canonicalized)
* @return set of bases to append to "node" to build the existing neighbors
*/
std::set<dna::base> UnitigBuilder::nextNeighbors(CQF<KeyObject> &cqf, dna::kmer& node) {
    std::set<dna::base> result;
    for (const auto b : dna::bases) {
        colorIdType eqid = 0;
        if (exists(cqf, node << b, eqid)) {
            result.insert(b);
        }
    }
    return result;
}

/**
 * search for a kmer in cqf and return the correct colorId if found
 * which is cqf count value - 1
 * @param cqf
 * @param e : search canonical kmer
 * @param eqid : reference to eqid that'll be set
 * @return true if eqid is found
 */
bool UnitigBuilder::exists(CQF<KeyObject> &cqf, dna::kmer e, colorIdType &eqid) {
    dna::canonical_kmer ck(e);
    KeyObject key(ck.val, 0, 0);
    auto eqidtmp = cqf.query(key, QF_NO_LOCK /*QF_KEY_IS_HASH | QF_NO_LOCK*/);
    if (eqidtmp) {
        eqid = eqidtmp - 1;
        return true;
    }
    return false;
}

void UnitigBuilder::buildUnitigs(CQF<KeyObject> &cqf) {
    logger->info("cqf.numslots: {}", cqf.numslots());
    logger->info("cqf.dist_elts: {}", cqf.dist_elts() );
    std::vector<bool> visited(cqf.numslots(), false);
    CQF<KeyObject>::Iterator it(cqf.begin());
    std::string filename(prefix + "/" + mantis::UNITIG_FASTA_FILE);
    std::ofstream out(filename, std::ios::out);
    uint64_t cntr{1}, visitedCnt{0};
    while (!it.done()) {
        KeyObject keyobj = *it;
        KeyObject key(keyobj.key, 0, 0);
        int64_t eqidx = cqf.get_unique_index(key, 0);
        if (!visited[eqidx]) { // for each unvisited kmer
            dna::kmer orig_node(static_cast<int>(k), keyobj.key);
            std::string seq = std::string(orig_node); // the unitig sequence
            dna::kmer curr_node = orig_node;
            visited[eqidx] = true;
            visitedCnt++;
            // extend it to the right as long as it has 1 and only 1 neighbor which is not visited yet
            auto neigbors = nextNeighbors(cqf, curr_node);
            while (neigbors.size() == 1) {
                dna::base n = (*neigbors.begin());
                curr_node = curr_node << n;
                dna::canonical_kmer cval(curr_node);
                eqidx = cqf.get_unique_index(KeyObject(cval.val, 0, 0), 0);
                if (visited[eqidx]) break;
                seq += dna::base_to_char.at(n);
                visited[eqidx] = true;
                visitedCnt++;

                neigbors = nextNeighbors(cqf, curr_node);
            }
            // extend it to the left with the same condition as right
            curr_node = orig_node;
            neigbors = prevNeighbors(cqf, curr_node);
            while (neigbors.size() == 1) {
                dna::base n = (*neigbors.begin());
                curr_node = n >> curr_node;
                dna::canonical_kmer cval(curr_node);
                eqidx = cqf.get_unique_index(KeyObject(cval.val, 0, 0), 0);
                if (visited[eqidx]) break;
                seq = dna::base_to_char.at(n) + seq;
                visited[eqidx] = true;
                visitedCnt++;
                neigbors = prevNeighbors(cqf, curr_node);
            }
            // write down the unitig with id u<cntr> and sequence <seq> in a fasta format
            out << ">u" << cntr << "\n";
            out << seq << "\n";
            cntr++;
        }
        ++it;
    }
    logger->info("total number of kmers (visited): {}", visitedCnt);
}

int main(int argc, char* argv[]) {
    auto console = spdlog::stdout_color_mt("mantis_console");
    if (argc < 2) {
        console->info("Usecase: buildUnitig <parent_dir of the cqf_file:{}> <num_of_threads (default:4)>", mantis::CQF_FILE);
        std::exit(1);
    }
    uint32_t numThreads = 4;
    if (argc == 3) {
        numThreads = static_cast<uint32_t >(std::stoul(argv[2]));
    }
    UnitigBuilder unitigBuilder(argv[1], console, numThreads);
}