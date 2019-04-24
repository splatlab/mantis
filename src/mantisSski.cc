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
    console->info("cnts: {}, {}, {}", cnt1, cnt2, cnt3);
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
//        std::cerr << "last one\n";
//        std::cerr << curBucketSize << " " << prevBucketstotalSize+curBucketSize << "\n";
        contigSeq.set_int((prevBucketstotalSize+curBucketSize)*2, 0, sizes.slicePrefixSize);
    }
    parser.stop();
    console->info("filled contig sequence vector and start idx vector.");
    console->info("cnts: {}, {}, {}", cnt1, cnt2, cnt3);
    // store the 2bit-encoded references
    sdsl::store_to_file(contigSeq, outdir + "/seq.bin");
    //sdsl::store_to_file(contigStartIdx, outdir + "/offset.bin");
}

/**
 * Builds the minimum perfect hash function
 * required information that have to be set in previous steps:
 * number of kmers (nkeys) -- contig sequence (contigSeq) and start index vector (contigStartIdx) -- kmer length (k)
 * k is limited to 2^8-1
 * @param numThreads
 */
void MantisSski::buildMPHF(uint32_t numThreads) {
//    std::cerr << "\n\nbuildMPHF\n\n\n";
//    std::cerr << " storage size: " << contigSeq.size() << "\n";
    ContigKmerIterator kb(&contigSeq, sizes, static_cast<uint8_t>(k), 0);
    ContigKmerIterator ke(&contigSeq, sizes, static_cast<uint8_t>(k), contigSeq.size() - k + 1);
    /*std::unordered_set<uint64_t> kmers;
    uint64_t kmrcntr=0;
    while (kb != ke) {
        auto t = *kb;
        if (kmers.find(t) != kmers.end()) {
            std::cerr << kmrcntr << " found it:\n";
            std::cerr << t << "\n";
            CanonicalKmer kk;
//            kk.fromNum(*kb);
//            std::cerr << kk.to_str() << "\n";
//            if (kmrcntr > 102480)
//            break;
        }
        kmrcntr++;
        kmers.insert(t);
        kb++;
    }
    std::exit(3);*/
    auto keyIt = boomphf::range(kb, ke);
    //std::cerr << " initialized the keyIt\n";
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
        auto contigLength = contigIdx+1!=contigCnt?
                contigStartIdx[contigIdx+1]-contigStartIdx[contigIdx]:
                contigSeq.size()-contigStartIdx[contigIdx];
        auto totalKmersIshouldSee = (contigLength - k + 1);
        contigIdx++ ;

        auto idx = bphf->lookup(*kb);
        // 0 : isCanonical, 1 otherwise
        // In case of stop char, use the reverse
        //We're at the beginning of the contig
        prefixArr[idx] = kb.isCanonical()? 1 : 0; // stop char
        kb++;
        cntr++;
        for (size_t j = 1; j < totalKmersIshouldSee; ++kb, ++j, ++cntr) {
            idx = bphf->lookup(*kb);
            int c = contigSeq[kb.pos()-1];
            if (kb.isCanonical()) {
                prefixArr[idx] = c << 1 & 0x6; // prefix char + 0 (for it is canonical)
            }
            else {
                c = combinelib::kmers::complement(c);
                prefixArr[idx] = c << 1 & 0x7; // prefix char + 1
            }
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