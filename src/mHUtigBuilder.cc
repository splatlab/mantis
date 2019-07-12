//
// Created by Fatemeh Almodaresi on 2019-05-10.
//

#include "mHUtigBuilder.h"

MHUtigBuilder::MHUtigBuilder(std::string prefixIn, std::shared_ptr<spdlog::logger> loggerIn, uint32_t numThreads) :
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
std::set<dna::base> MHUtigBuilder::prevNeighbors(CQF<KeyObject> &cqf, dna::kmer& node) {
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
std::set<dna::base> MHUtigBuilder::nextNeighbors(CQF<KeyObject> &cqf, dna::kmer& node) {
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
bool MHUtigBuilder::exists(CQF<KeyObject> &cqf, dna::kmer e, colorIdType &eqid) {
    dna::canonical_kmer ck(e);
    KeyObject key(ck.val, 0, 0);
    auto eqidtmp = cqf.query(key, QF_NO_LOCK /*QF_KEY_IS_HASH | QF_NO_LOCK*/);
    if (eqidtmp) {
        eqid = eqidtmp - 1;
        return true;
    }
    return false;
}

void MHUtigBuilder::buildUnitigs(CQF<KeyObject> &cqf) {
    logger->info("cqf.numslots: {}", cqf.numslots());
    logger->info("cqf.dist_elts: {}", cqf.dist_elts() );
    std::vector<bool> visited(cqf.numslots(), false);
    CQF<KeyObject>::Iterator it(cqf.begin());
    std::string filename(prefix + "/" + mantis::MHUTIG_FASTA_FILE);
    logger->info("output file name: {}", filename);
    std::ofstream out(filename, std::ios::out);
    uint64_t cntr{1}, visitedCnt{0};
    while (!it.done()) {
        KeyObject keyobj = *it;
        KeyObject key(keyobj.key, 0, 0);
        int64_t keyidx = cqf.get_unique_index(key, 0);
//        std::cerr << "next round: " << keyidx << "\n";
        if (!visited[keyidx]) { // for each unvisited kmer
            dna::kmer orig_node(static_cast<int>(k), keyobj.key);
            std::string seq = std::string(orig_node); // the unitig sequence
            dna::kmer curr_node = orig_node;
            visited[keyidx] = true;
            visitedCnt++;

            while (true) {
                auto neigbors = nextNeighbors(cqf, curr_node);
                std::vector<std::pair<uint64_t, dna::base>> charList;
                charList.reserve(8);
                for (auto n : neigbors) {
                    auto tmp_node = curr_node << n;
                    dna::canonical_kmer cval(tmp_node);
                    keyidx = cqf.get_unique_index(KeyObject(cval.val, 0, 0), 0);
                    if (!visited[keyidx]) charList.emplace_back(keyidx, n);
                }
                if (charList.empty()) break;
//                std::cerr << "n" << charList.size() << "\n";
                // randomly select a base
                auto k = std::rand() % charList.size();
//                std::cerr << "k" << k << "\n";
                dna::base n = charList[k].second;
                seq += dna::base_to_char.at(n);
                visited[charList[k].first] = true;
                visitedCnt++;
                curr_node = curr_node << n;
            }
//            std::cerr << "\nprevious\n";
            curr_node = orig_node;
            while (true) {
//                std::cerr << "Start\n";
                auto neigbors = prevNeighbors(cqf, curr_node);
                std::vector<std::pair<uint64_t, dna::base>> charList;
                charList.reserve(8);
                for (auto nItr = neigbors.begin(); nItr != neigbors.end(); nItr++) {
                    dna::base n = (*nItr);
                    auto tmp_node = n >> curr_node;
                    dna::canonical_kmer cval(tmp_node);
                    keyidx = cqf.get_unique_index(KeyObject(cval.val, 0, 0), 0);
                    if (!visited[keyidx]) charList.emplace_back(keyidx, n);
                }
                if (charList.empty()) break;
//                std::cerr << "n" << charList.size() << "\n";
                // randomly select a base
                auto k = std::rand() % charList.size();
//                std::cerr << "k" << k << "\n";
                dna::base n = charList[k].second;
//                std::cerr << seq << " " << dna::base_to_char.at(n) << "\n";
                seq = dna::base_to_char.at(n) + seq;
                visited[charList[k].first] = true;
                visitedCnt++;
                curr_node = n >> curr_node;
            }
//            std::cerr << "1\n";
            // write down the unitig with id u<cntr> and sequence <seq> in a fasta format
            out << ">u" << cntr << "\n";
            out << seq << "\n";
            cntr++;
//            std::cerr << "2\n";
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
    MHUtigBuilder unitigBuilder(argv[1], console, numThreads);
}