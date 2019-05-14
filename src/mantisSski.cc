//
// Created by Fatemeh Almodaresi on 2019-03-07.
//

#include <ProgOpts.h>
#include "mantisSski.h"
#include "FastxParser.hpp"
#include <unordered_set>

/**
 *
 * @param seqVec contigs sequence vector
 * @param offset the position contig sequence should be put in seqVec
 * @param str the contig sequence
 */
void MantisSski::encodeSeq(sdsl::int_vector<2> &seqVec, size_t offset, stx::string_view str) {
//    std::cerr << "encodeSeq: " << offset << ", sv:" << str.to_string() << "\n";
    for (size_t i = 0; i < str.length(); ++i) {
        auto c = combinelib::kmers::codeForChar(str[i]);
//        std::cerr << c;
        seqVec[offset + i] = c;
    }
//    std::cerr << "\n";
}


/**
 * Parses the unitig file twice.
 * Calculates the overal length of the contig vector and number of kmers in the first round
 * Converts contigs to 2-bit encoded sequences and puts them in the vector
 * @param numThreads
 * @param cfile contig file
 */
void MantisSski::buildUnitigVec(uint32_t numThreads, std::string cfile) {
    // go over all the contig files
    console->info("Reading the contig file ...");
    std::vector<std::string> ref_files = {cfile};
    //FIXME how to use one parser to parse the file twice
    fastx_parser::FastxParser<fastx_parser::ReadSeq> parser1(ref_files, numThreads, 1);
    parser1.start();
    auto rg = parser1.getReadGroup();
    // read the reference sequences and encode them into the refseq int_vector
    uint64_t bucketSizeInInt{sizes.bucketSize/2}, slicePrefixSizeInInt{sizes.slicePrefixSize/2};
    uint64_t contigAccumLength{bucketSizeInInt}, remaining{0}, quotient{0}, curBucketSize{0}, added2bits{0};
    uint64_t cnt1{0}, cnt2{0}, cnt3{0};
    while (parser1.refill(rg)) {
        for (auto &rp : rg) {
            cnt1++;
            std::string seq = rp.seq;
            // this is the ridiculous condition for taking care of a corner case in bcalm
            if (seq.size() > k and seq.substr(0, k) == seq.substr(seq.size()-k, k)) {
                seq = seq.substr(0, seq.size()-1);
            }
            uint64_t curSize = seq.size();
//            if (rp.seq.size() > k and isPalyndrome(rp.seq))
//                std::cout << rp.seq << "\n";
            do {
                cnt2++;
                remaining = curSize % (sizes.unitigSliceSize + 1);
                quotient = curSize / (sizes.unitigSliceSize + 1);
                curSize = curSize - sizes.unitigSliceSize + k-1; // update curSize for next iteration
                if (quotient) {
                    added2bits = sizes.unitigSliceSize;
                } else {
                    added2bits = remaining;
                }
                if (curBucketSize + added2bits + slicePrefixSizeInInt > bucketSizeInInt) {
                    cnt3++;
                    contigAccumLength += bucketSizeInInt;
                    curBucketSize = 0;
                }
                curBucketSize += added2bits + slicePrefixSizeInInt; // for next round. If last round, then we don't need it
            } while (quotient);

            nkeys += seq.size() - k + 1;
            contigCnt++;
        }
    }
    parser1.stop();
    console->info("# of contigs: {}", contigCnt);
    console->info("# of kmers: {}", nkeys);
    console->info("contigSeq Length: {}", contigAccumLength);
    console->info("# of contigs: {}, # of 256B slices: {}, # of 1k buckets: {}", cnt1, cnt2, cnt3);
    // second round over the file
    contigSeq = sdsl::int_vector<2>(contigAccumLength, 0);
    //contigStartIdx = sdsl::int_vector<>(contigCnt, 0, std::log2(contigAccumLength) + 1);
    fastx_parser::FastxParser<fastx_parser::ReadSeq> parser(ref_files, numThreads, 1);
    parser.start();
    rg = parser.getReadGroup();
    // read the reference sequences and encode them into the refseq int_vector
    uint64_t prevBucketstotalSize{0}, bucketContigCntr{0}, contigCntr{0}, startIdx;
    curBucketSize=0;remaining=0;quotient=0;
    cnt1=0;cnt2=0;cnt3=0;
    while (parser.refill(rg)) {
        for (auto &rp : rg) {
            cnt1++;
            startIdx=0;
            std::string seq = rp.seq;

            if (seq.size() > k and seq.substr(0, k) == seq.substr(seq.size()-k, k)) {
                seq = seq.substr(0, seq.size()-1);
            }
            stx::string_view seqv(seq);
            uint64_t curSize = seq.size();

            do {
                cnt2++;
                remaining = curSize % (sizes.unitigSliceSize + 1);
                quotient = curSize / (sizes.unitigSliceSize + 1);
                curSize = curSize - sizes.unitigSliceSize + k-1;
                if (quotient) {
                    added2bits = sizes.unitigSliceSize;
                } else {
                    added2bits = remaining;
                }
                stx::string_view subsv = seqv.substr(startIdx, added2bits);
                startIdx += (added2bits-k+1);
                if (curBucketSize + added2bits + slicePrefixSizeInInt > bucketSizeInInt) {
                    cnt3++;
                    // # of bits remaining < slicePrefixSize
                    // or
                    // a slice of size zero are the two signatures for the end of the bucket
                    if (bucketSizeInInt - curBucketSize >= slicePrefixSizeInInt) {
                        contigSeq.set_int((prevBucketstotalSize+curBucketSize)*2, 0, sizes.slicePrefixSize);
                    }
                    bucketContigCntr = 0;
                    prevBucketstotalSize += bucketSizeInInt;
                    curBucketSize = 0;
                }
                contigSeq.set_int((prevBucketstotalSize+curBucketSize)*2, subsv.size(), sizes.slicePrefixSize);
                curBucketSize+=slicePrefixSizeInInt;
                encodeSeq(contigSeq, prevBucketstotalSize+curBucketSize, subsv);
                curBucketSize += added2bits;
            } while (quotient);
            contigCntr++;
            bucketContigCntr++;
        }
    }
    if (bucketSizeInInt - curBucketSize >= slicePrefixSizeInInt) {
        contigSeq.set_int((prevBucketstotalSize+curBucketSize)*2, 0, sizes.slicePrefixSize);
    }
    parser.stop();
    console->info("filled contig sequence vector and start idx vector.");
    console->info("cnts: {}, {}, {}", cnt1, cnt2, cnt3);
    // store the 2bit-encoded references
    sdsl::store_to_file(contigSeq, outdir + "/seq.bin");
}

/**
 * Builds the minimum perfect hash function
 * required information that have to be set in previous steps:
 * number of kmers (nkeys) -- contig sequence (contigSeq) and start index vector (contigStartIdx) -- kmer length (k)
 * k is limited to 2^8-1
 * @param numThreads
 */
void MantisSski::buildMPHF(uint32_t numThreads) {
    ContigKmerIterator kb(&contigSeq, sizes, static_cast<uint8_t>(k), 0);
    ContigKmerIterator ke(&contigSeq, sizes, static_cast<uint8_t>(k), contigSeq.size() - k + 1);
    auto keyIt = boomphf::range(kb, ke);
    bphf =
            new boophf_t(nkeys, keyIt, numThreads, 3.5); // keys.size(), keys, 16);
    console->info("mphf size = {}", (bphf->totalBitSize() / 8) / std::pow(2, 20));
    // store to disk
    std::ofstream hstream(outdir + "/mphf.bin");
    bphf->save(hstream);
    hstream.close();
    console->info("stored mphf.");
}

/**
 * Builds prefixArr, a 3-bit array that contains either the char to prepend to get to the previous contig
 * or stop char if we are at the beginning
 * Since a kmer is not necessarily canonicalized in a unitig, we need one extra bit to store if kmer is canonicalized
 */
void MantisSski::buildPrefixArr() {
    console->info("building the prefix array ...");

    prefixArr = sdsl::int_vector<3>(nkeys, 0);
    ContigKmerIterator kb(&contigSeq, sizes, static_cast<uint8_t>(k), 0);
    ContigKmerIterator ke(&contigSeq, sizes, static_cast<uint8_t>(k), contigSeq.size() - k + 1);

    uint64_t contigIdx{0};
    // For every valid k-mer (i.e. every contig)
    uint64_t cntr{0}, prevPrcnt{0}, curPrcnt{0};
    while(kb != ke){
        auto totalKmersIshouldSee = kb.getSliceSize()-k+1;
        contigIdx++ ;

        auto idx = bphf->lookup(*kb);
        // 0 : isCanonical, 1 otherwise
        // In case of stop char, use the reverse
        //We're at the beginning of the contig
        if (totalKmersIshouldSee == 1) {
            prefixArr[idx] = 0;
        } else {
            int c = contigSeq[kb.pos()+k];
            c = combinelib::kmers::complement(c);
            // for the first kmer in the contig, the canonical bit should show the reverse
            // so that it creates the next kmer
            prefixArr[idx] = c | (!kb.isCanonical() << 2); // first 2 bits are for the nucleotide
        }
        kb++;
        cntr++;
        for (size_t j = 1; j < totalKmersIshouldSee; ++kb, ++j, ++cntr) {
            idx = bphf->lookup(*kb);
            int c = contigSeq[kb.pos()-1];
            prefixArr[idx] = c | (kb.isCanonical() << 2); // first 2 bits are for the nucleotide
            /*if (kb.isCanonical()) {
                std::cerr << (uint16_t)kb.isCanonical() << " "
                << (uint16_t)c << " "
                << (uint16_t)(kb.isCanonical() << 2) << " " << (uint16_t)(c | (kb.isCanonical() << 2))<<"\n";
                std::exit(3);
            }*/
        }
        curPrcnt = (cntr*100)/nkeys;
        if (curPrcnt > prevPrcnt) {
            std::cerr << "\r " << curPrcnt << "% completed ...";
            prevPrcnt = curPrcnt;
        }
    }
    std::cerr << "\n";

    // store prefixArr
    console->info("storing the prefixArray ...");
    sdsl::store_to_file(prefixArr, outdir + "/prefix.bin");
    console->info("stored prefixArray.");
}

bool MantisSski::queryKmer(CanonicalKmer kmer) {
    /*auto numBuckets = contigSeq.size()/(sizes.bucketSize/2);
    for (auto mid=0; mid < numBuckets ; mid++) {
        uint64_t midWrd = contigSeq.get_int(sizes.bucketSize * mid + sizes.slicePrefixSize, 2 * k);
        CanonicalKmer k;
        k.fromNum(midWrd);
        std::cerr << "bucket " << mid << ": " << midWrd << " " << k.to_str() << "\n";
    }
    std::exit(3);*/
    Cases cases = Cases::last2kmers;
    // conflicted: while condition --> stop going through mphf as soon as you reach a conflict (which happens in multiple cases)
    // first: just a flag for the first round of MPHF+sequenceVecs
    bool first = true;
    // firstPrefix: prefix of the first k-mer --> interesting if 000 in case of conflict
    int firstPrefix{1}, cntr{0};
    CanonicalKmer originalKmer{kmer}, prevPrevKmer, prevKmer;
    prevKmer.fromNum(0);
    prevPrevKmer.fromNum(0);
    std::cerr << "prefix.size(): " << prefixArr.size() << "\n";
    while (true) {
        std::cerr << kmer.to_str() << "\n";
        auto idx = bphf->lookup(kmer.getCanonicalWord());
        std::cerr << "idx:" << idx << " ";
        // if the returning number by MPHF is out of range or we've exhausted the size of a slice
        //todo check corner case : sliceSize - k or sliceSize - k + 1
        if (idx > prefixArr.size() or cntr > sizes.slicePrefixSize - k) {
            std::cerr << "got to a conflicting point. idx > prefixArr.size()? "
            << (bool)(idx > prefixArr.size()) << "\n";
            if (!firstPrefix) {
                cases =  Cases::firstKmer;
                break;
            }
            else {
                cases = Cases::notFound;
                break;
            }
        }
        // attach prefix
        uint8_t val = prefixArr[idx];
        std::cerr << idx << " prefix val: " <<  (uint16_t)val << "\n";
        int prefix = val & 0x3;
        if (first) {
            firstPrefix = val;
            first = false;
        }
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
            cases = Cases::last2kmers;
            break;
        }
        cntr++;
    }
    if (cases == Cases::notFound) {
        std::cerr << "notFound\n";
        return false;
    } else if (cases == Cases::last2kmers) {
        std::cerr << "last2kmers\n";
        auto idx = binarySearch(kmer);
        if (searchBucket(kmer, idx)) return true;
        kmer.swap();
        idx = binarySearch(kmer);
        if (searchBucket(kmer, idx)) return true;

        idx = binarySearch(prevKmer);
        if (searchBucket(prevKmer, idx)) return true;
        prevKmer.swap();
        idx = binarySearch(prevKmer);
        if (searchBucket(prevKmer, idx)) return true;
        return searchBucket(prevKmer, idx);
    } else {
        std::cerr << "firstkmer\n";
        auto idx = binarySearch(originalKmer);
        if (searchBucket(originalKmer, idx)) return true;
        originalKmer.swap();
        idx = binarySearch(originalKmer);
        return searchBucket(originalKmer, idx);
    }
    /*switch (cases) {
        case Cases::notFound:
            return false;
        case Cases::last2kmers:
            uint64_t idx = binarySearch(kmer);
            if (searchBucket(kmer, idx)) return true;
            idx = binarySearch(prevKmer);
            return searchBucket(prevKmer, idx);
        case Cases::firstKmer:
            uint64_t idx = binarySearch(originalKmer);
            return searchBucket(originalKmer, idx);
    }*/
}

uint64_t MantisSski::binarySearch(CanonicalKmer kmer) {
    std::cerr << "Binary Search\n";
    std::cerr << "2bSearched: " << kmer.to_str() << "\n";
    uint64_t wrd = kmer.fwWord();
        auto numBuckets = contigSeq.size()/(sizes.bucketSize/2);
        uint64_t low{0}, high{numBuckets - 1};
        std::cerr << numBuckets << "\n";
        // The algorithm is a little bit different
        // we don't set high and low index to mid+1 and mid-1 respectively.
        // Because we're not searching for the exact kmer but a range (The correct bucket)
        while (1 < high - low) {
            auto mid = (low + high) / 2;
            // go to the "mid" bucket. fetch a word from after the sliceSize prefix
            uint64_t midWrd = contigSeq.get_int(sizes.bucketSize*mid+sizes.slicePrefixSize, 2*k);
            CanonicalKmer k;
            k.fromNum(midWrd);
            std::cerr << "bucket " << mid << ": " << k.to_str() << " ";
            if (midWrd > wrd) {
                std::cerr << "midWrd > wrd " << midWrd << ">" << wrd << "\n";
                /*for (int i = 0; i < 47; i++) {
                    std::cerr << ((midWrd >> i) & 0x01)?"1":"0";
                }
                std::cerr << "\n";*/
                high = mid;
            } else if (midWrd < wrd) {
                std::cerr << "midWrd < wrd " << midWrd << "<" << wrd << "\n";
                low = mid;
            } else {
                std::cerr << "midWrd==wrd\n";
                return mid;
            }
        }
        return low;
}

bool MantisSski::searchBucket(CanonicalKmer kmer, uint64_t idx) {
    uint64_t bstart{sizes.bucketSize*idx}, bend{(sizes.bucketSize+1)*idx};
    std::cerr << "Bucket Search:\n";
    std::cerr << "2bSearched: " << kmer.to_str() << "\n";
    uint64_t cntr = 0;
    while (bend - bstart > sizes.slicePrefixSize) {
        auto sliceLength = contigSeq.get_int(bstart, sizes.slicePrefixSize);
        bstart+= sizes.slicePrefixSize;
        if (!sliceLength) break;
        uint64_t wrd = contigSeq.get_int(bstart, 2*k);
        CanonicalKmer searched;
        searched.fromNum(wrd);
        std::cerr << cntr++ << " " << sliceLength << " " << searched.to_str() <<  " " << searched.fwWord() <<"\n";
        if (searched.getCanonicalWord() == kmer.getCanonicalWord()) return true;
        bstart+=(2*sliceLength);
    }
    return false;
}

int build_sski_main ( BuildOpts& opt )
{
    MantisSski mantisSski(opt.k, opt.out, opt.console.get());
    // these three should be called in order
    mantisSski.buildUnitigVec(static_cast<uint32_t>(opt.numthreads), opt.inlist);
//    mantisSski.testBooPHF(opt.numthreads);
    mantisSski.buildMPHF(static_cast<uint32_t>(opt.numthreads));
    mantisSski.buildPrefixArr();
    return 0;
}

/*
 * ===  FUNCTION  =============================================================
 *         Name:  main
 *  Description:
 * ============================================================================
 */
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
    MantisSski sski(queryK, opt.prefix, logger, true);
    logger->info("Done loading sski.");

    //std::ofstream opfile(opt.output);
    std::ifstream ipfile(opt.query_file);
    std::string read;
    uint64_t numOfQueries{0};
    CanonicalKmer kmer;
    while (ipfile >> read) {
        std::cerr << read << " ";
        kmer.fromStr(read);
        std::cerr << sski.queryKmer(kmer) << " -- > Found the kmer " << kmer.to_str() << "\n";
        numOfQueries++;
    }
    /*CLI::AutoTimer timer{"query time ", CLI::Timer::Big};
    if (opt.process_in_bulk) {
        while (ipfile >> read) {
            mstQuery.parseKmers(read, indexK);
            numOfQueries++;
        }
        mstQuery.findSamples(cqf, cache_lru, &rs, queryStats);
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
    }*/
    return EXIT_SUCCESS;
}
