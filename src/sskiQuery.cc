//
// Created by Fatemeh Almodaresi on 2019-07-12.
//

#include "sskiQuery.h"
#include <ProgOpts.h>
#include <CLI/Timer.hpp>

#define BITMASK(nbits) ((nbits) == 64 ? 0xffffffffffffffff : (1ULL << (nbits)) - 1ULL)


bool SskiQuery::queryKmer(CanonicalKmer kmer) {

    // first: just a flag for the first round of MPHF+sequenceVecs
    bool first = true;
    // firstPrefix: prefix of the first k-mer --> interesting if 000 in case of conflict
    int cntr{0};
    queriedKmer = kmer;
    CanonicalKmer originalKmer{kmer}, prevPrevKmer, prevKmer;
    prevKmer.fromNum(0);
    prevPrevKmer.fromNum(0);

    stepsFw = 0;
    while (true) {
        auto idx = bphf->lookup(kmer.getCanonicalWord());
        // if the returning number by MPHF is out of range or we've exhausted the size of a slice
        if (idx > prefixArr.size() or cntr > sizes.unitigSliceSize - k) {
//            std::cerr << "mphf couldn't find it\n";
            return false;
        }
        // fetch prefix
        uint8_t val = prefixArr[idx];
        int prefix = val & 0x3;
        // If it is the first kmer to search and the prefix value is 0, It can be a special case for uniq-kmer contigs
        // Go search for such a contig
        if (first and val == 0) {
            auto firstIdx = binarySearch(originalKmer);
            if (searchBucket(originalKmer, firstIdx)) return true;
            originalKmer.swap();
            firstIdx = binarySearch(originalKmer);
            if (searchBucket(originalKmer, firstIdx)) return true;
        }
        first = false;
        prevPrevKmer = prevKmer;
        prevKmer = kmer;
        if (val & 0x4) { // go for canonical kmer in the seq vec
            if (kmer.isFwCanonical()) // if the kmer you're searching is also canonical
                kmer.shiftBw(prefix);
            else
                kmer.shiftFw(combinelib::kmers::complement(prefix));
        } else { // go for the non-canonical version of kmer
            if (kmer.isFwCanonical()) // if the kmer you're searching is also non-canonical
                kmer.shiftFw(combinelib::kmers::complement(prefix));
            else
                kmer.shiftBw(prefix);
        }
        if (kmer == prevPrevKmer) {
            break;
        }
//        std::cerr << prevPrevKmer.to_str() << " " << prevKmer.to_str() << " " << kmer.to_str() << "\n";
        cntr++;
    }
    stepsFw= static_cast<uint64_t>(cntr);

    // You reached to a case that the last two kmers point at each other
    // Search for both of them
    auto idx = binarySearch(prevKmer); // prevKmer fw
    if (searchBucket(prevKmer, idx)) return true;
    prevKmer.swap();
    idx = binarySearch(prevKmer); // prevKmer rc
    if (searchBucket(prevKmer, idx)) return true;

    idx = binarySearch(kmer); // kmer fw
    if (searchBucket(kmer, idx)) return true;
    kmer.swap();
    idx = binarySearch(kmer); // kmer rc
    if (searchBucket(kmer, idx)) return true;
}

uint64_t SskiQuery::binarySearch(CanonicalKmer kmer) {
    uint64_t wrd = kmer.fwWord();
    auto numBuckets = contigSeq.size() / (sizes.bucketSize / 2);
    uint64_t low{0}, high{numBuckets};
//    std::cerr << "binary search " << kmer.to_str() << "\n";
    // The algorithm is a little bit different
    // we don't set high and low index to mid+1 and mid-1 respectively.
    // Because we're not searching for the exact kmer but a range (The correct bucket)
    while (1 < high - low) {
        auto mid = (low + high) / 2;
        // go to the "mid" bucket. fetch a word from after the sliceSize prefix
        uint64_t midWrd = contigSeq.get_int(sizes.bucketSize * mid + sizes.slicePrefixSize, 2 * k);
        CanonicalKmer k;
        k.fromNum(midWrd);
//        std::cerr << mid << " " << k.to_str() << "\n";
        if (midWrd > wrd) {
            high = mid;
        } else if (midWrd < wrd) {
            low = mid;
        } else {
//            std::cerr << "return mid: " << mid << "\n";
            return mid;
        }
    }
//    std::cerr << "return low: " << low << "\n";
    return low;
}

uint64_t SskiQuery::searchBucket(CanonicalKmer kmer, uint64_t idx) {
    uint64_t bstart{sizes.bucketSize * idx}, bend{(sizes.bucketSize) * (idx + 1)};
//    std::cerr << "search bucket\n";
    while (bend > bstart and bend - bstart > sizes.slicePrefixSize) {
        auto sliceLength = contigSeq.get_int(bstart, sizes.slicePrefixSize);
//        std::cerr << bstart << " " << sliceLength << " ";
        bstart += sizes.slicePrefixSize;
        if (!sliceLength) break;
        uint64_t wrd = contigSeq.get_int(bstart, 2 * k);
        CanonicalKmer searched;
        searched.fromNum(wrd);
//        std::cerr << searched.to_str() << "\n";
        // bstart will never be 0 because we return the index after the prefix ends
        if (searched.getCanonicalWord() == kmer.getCanonicalWord()) {
//            std::cerr << "steps: " << stepsFw << "\n";
            if (stepsFw+k-1<=sliceLength) {
                auto origWrd = contigSeq.get_int(bstart + (stepsFw - 1) * 2, 2 * k);
                CanonicalKmer kk;
                kk.fromNum(origWrd);
//                std::cerr << "origWrd: " << kk.to_str() << "\n";
                if (origWrd == queriedKmer.fwWord() or origWrd == queriedKmer.rcWord())
                    return bstart + (stepsFw - 1) * 2;
            }
            if (stepsFw+k<=sliceLength) {
                auto origWrd = contigSeq.get_int(bstart + stepsFw * 2, 2 * k);
                CanonicalKmer kk;
                kk.fromNum(origWrd);
//                std::cerr << "origWrd: " << kk.to_str() << "\n";
                if (origWrd == queriedKmer.fwWord() or origWrd == queriedKmer.rcWord())
                    return bstart + stepsFw * 2;
            }
            // it was a fp
            return 0;
            return bstart;
        }
        bstart += (2 * sliceLength);
    }
    return 0; // 0 means not found
}


void SskiQuery::parseKmers(std::string read, uint64_t kmer_size) {
    //CLI::AutoTimer timer{"First round going over the file ", CLI::Timer::Big};
    bool done = false;
    while (!done and read.length() >= kmer_size) {
        uint64_t first = 0;
        uint64_t first_rev = 0;
        uint64_t item = 0;
        bool allNuclValid = true;
        for (uint32_t i = 0; i < kmer_size; i++) { //First kmer
            int curr = kmers::codeForChar(read[i]);
            if (curr == -1) { // 'N' is encountered
                if (i + 1 < read.length())
                    read = read.substr(i + 1);
                else
                    done = true; // you've reached the end of a sequence (last nucl. is N)
                allNuclValid = false;
                break;
            }
            first = first | curr;
            first = first << 2;
        }
        if (allNuclValid) {
            first = first >> 2;
            first_rev = static_cast<uint64_t >(kmers::word_reverse_complement(first, static_cast<uint16_t>(kmer_size)));

            //cout << "kmer: "; cout << int_to_str(first);
            //cout << " reverse-comp: "; cout << int_to_str(first_rev) << endl;

            if (first < first_rev)
                item = first;
            else
                item = first_rev;
            kmer2cidMap[item] = std::numeric_limits<uint64_t>::max();

            uint64_t next = (first << 2) & BITMASK(2 * kmer_size);
            uint64_t next_rev = first_rev >> 2;
            uint64_t i = 0;
            for (i = kmer_size; i < read.length(); i++) { //next kmers
                //cout << "K: " << read.substr(i-K+1,K) << endl;
                int curr = kmers::codeForChar(read[i]);
                if (curr == -1) { // 'N' is encountered
                    if (i + 1 < read.length())
                        read = read.substr(i + 1);
                    else
                        done = true;
                    break;
                }
                next |= curr;
                auto tmp = static_cast<uint64_t>(kmers::revCodes[curr]);
                tmp <<= (kmer_size * 2 - 2);
                next_rev = next_rev | tmp;
                if (next < next_rev)
                    item = next;
                else
                    item = next_rev;
                kmer2cidMap[item] = std::numeric_limits<uint64_t>::max();
                next = (next << 2) & BITMASK(2 * kmer_size);
                next_rev = next_rev >> 2;
            }
            if (i == read.length()) done = true;
        }
    }
}

/*
bool MantisSski::collectStat() {
    std::ofstream overlaps("overlaps.out", "w");
    uint64_t cntr{0};
    uint64_t bstart{sizes.bucketSize * idx}, bend{(sizes.bucketSize) * (idx + 1)};
    auto sliceLength = contigSeq.get_int(bstart, sizes.slicePrefixSize);
    bstart += sizes.slicePrefixSize;
    uint64_t wrd = contigSeq.get_int(bstart, 2 * k);
    CanonicalKmer first;
    first.fromNum(wrd);
    auto prevSliceLength = sliceLength;
    uint16_t overlap{0};
    while (bend > bstart and bend - bstart > sizes.slicePrefixSize) {
        prevSliceLength = sliceLength;
        auto sliceLength = contigSeq.get_int(bstart, sizes.slicePrefixSize);
        if (!sliceLength) break;
        bstart += sizes.slicePrefixSize + (2 * sliceLength);
    }
    if (sliceLength) sliceLength = prevSliceLength;

    bstart -= sliceLength;
    uint64_t wrd = contigSeq.get_int(bstart, 2 * k);
    CanonicalKmer last;
    last.fromNum(wrd);

    ofstream << cntr << " " << overlap << "\n";
    return false;
}
*/


/*
 * ===  FUNCTION  =============================================================
 *         Name:  main
 *  Description:
 * ============================================================================
 */

#include "gqf/hashutil.h"
#include "gqf_cpp.h"

int sski_query_main(QueryOpts &opt) {
    uint32_t queryK = opt.k;
//    QueryStats queryStats;

    spdlog::logger *logger = opt.console.get();
    /*std::string dbg_file(opt.prefix + mantis::CQF_FILE);
    std::string sample_file(opt.prefix + mantis::SAMPLEID_FILE);

    std::vector<std::string> sampleNames = loadSampleFile(sample_file);
    queryStats.numSamples = sampleNames.size();
    logger->info("Number of experiments: {}", queryStats.numSamples);*/

    logger->info("Loading sski...");
    SskiQuery sski(queryK, opt.prefix, logger, true);
    logger->info("Done loading sski.");

//    std::ofstream opfile(opt.output);
    std::ifstream ipfile(opt.query_file);
    std::ofstream out(opt.output);
    std::string read;
    uint64_t numOfQueries{0};
    //CanonicalKmer kmer;

//    std::string cqf_file(opt.query_file);
//    CQF<KeyObject> cqf(cqf_file, CQF_FREAD);
//    logger->info("Done loading the cqf file.");
//    auto k = cqf.keybits() / 2;
//    auto it = cqf.begin();
    CLI::AutoTimer timer{"query time ", CLI::Timer::Big};
    while (ipfile >> read) {
//    while (!it.done()) {
//        KeyObject keyobj = *it;
//        dna::canonical_kmer kmertmp(k, keyobj.key);
        CanonicalKmer kmer;
//        kmer.fromStr(std::string(kmertmp));
        kmer.fromStr(read);
        /*if (!sski.queryKmer(kmer)) {
            std::cerr << "\n" << numOfQueries << " " << kmer.to_str() << "\n";
//            std::exit(3);
        }*/
        auto res = sski.queryKmer(kmer);
        if (res)
            out << kmer.to_str() << " f\n";
        else
            out << kmer.to_str() << " nf\n";
        numOfQueries++;
//        ++it;
        if (numOfQueries % 1000000 == 0)
            std::cerr << "\r" << numOfQueries;
    }
    std::cerr << "\nNum of queries: " << numOfQueries << "\n";
    out.close();
/*

    CLI::AutoTimer timer{"query time ", CLI::Timer::Big};
    if (opt.process_in_bulk) {
        while (ipfile >> read) {
            sski.parseKmers(read, k);
            numOfQueries++;
        }
        sski.findSamples(cqf, cache_lru, &rs, queryStats);
        ipfile.clear();
        ipfile.seekg(0, ios::beg);
        if (opt.use_json) {
            opfile << "[\n";
            while (ipfile >> read) {
                output_results_json(read, mstQuery, opfile, sampleNames, queryStats, numOfQueries);
            }
            opfile << "]\n";
        } else {
            while (ipfile >> read) {
                output_results(read, mstQuery, opfile, sampleNames, queryStats);
            }
        }
    } else {
        if (opt.use_json) {
            opfile << "[\n";
            while (ipfile >> read) {
                mstQuery.reset();
                mstQuery.parseKmers(read, indexK);
                mstQuery.findSamples(cqf, cache_lru, &rs, queryStats);
                if (mstQuery.indexK == mstQuery.queryK)
                    output_results_json(mstQuery, opfile, sampleNames, queryStats, numOfQueries);
                else
                    output_results_json(read, mstQuery, opfile, sampleNames, queryStats, numOfQueries);
                numOfQueries++;
            }
            opfile << "]\n";
        } else {
            while (ipfile >> read) {
                mstQuery.reset();
                mstQuery.parseKmers(read, indexK);
                mstQuery.findSamples(cqf, cache_lru, &rs, queryStats);
                if (mstQuery.indexK == mstQuery.queryK)
                    output_results(mstQuery, opfile, sampleNames, queryStats);
                else
                    output_results(read, mstQuery, opfile, sampleNames, queryStats);
                numOfQueries++;
            }
        }
    }
    opfile.close();
    logger->info("Writing done.");

    for (auto &kv : queryStats.numOcc) {
        std::cout << kv.first << '\t' << kv.second << '\n';
    }
*/
    return EXIT_SUCCESS;
}
