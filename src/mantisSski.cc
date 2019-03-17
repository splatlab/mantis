//
// Created by Fatemeh Almodaresi on 2019-03-07.
//

#include <ProgOpts.h>
#include "mantisSski.h"
#include "FastxParser.hpp"

/**
 *
 * @param seqVec contigs sequence vector
 * @param offset the position contig sequence should be put in seqVec
 * @param str the contig sequence
 */
void MantisSski::encodeSeq(sdsl::int_vector<2> &seqVec, size_t offset, stx::string_view str) {
//    std::cerr << offset << ", sv:" << str.to_string() << "\n";
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
    uint64_t contigAccumLength{0};
    while (parser1.refill(rg)) {
        for (auto &rp : rg) {
            contigAccumLength+=rp.seq.size();
            nkeys += rp.seq.size() - k + 1;
            contigCnt++;
        }
    }
    parser1.stop();
    console->info("# of contigs: {}", contigCnt);
    console->info("# of kmers: {}", nkeys);

    // second round over the file
    contigSeq = sdsl::int_vector<2>(contigAccumLength, 0);
    contigStartIdx = sdsl::int_vector<>(contigCnt, 0, std::log2(contigAccumLength) + 1);
    fastx_parser::FastxParser<fastx_parser::ReadSeq> parser(ref_files, numThreads, 1);
    parser.start();
    rg = parser.getReadGroup();
    // read the reference sequences and encode them into the refseq int_vector
    uint64_t offset{0}, contigCntr{0};
    while (parser.refill(rg)) {
        for (auto &rp : rg) {
            stx::string_view seqv(rp.seq);
            encodeSeq(contigSeq, offset, seqv);
            contigStartIdx[contigCntr] = offset;
            //std::cerr << contigCntr << ":" << offset << "\n";
            contigCntr++;
            offset += rp.seq.size();
        }
    }
    parser.stop();
    console->info("filled contig sequence vector and start idx vector.");
    // store the 2bit-encoded references
    sdsl::store_to_file(contigSeq, outdir + "/seq.bin");
    sdsl::store_to_file(contigStartIdx, outdir + "/offset.bin");
}

/**
 * Builds the minimum perfect hash function
 * required information that have to be set in previous steps:
 * number of kmers (nkeys) -- contig sequence (contigSeq) and start index vector (contigStartIdx) -- kmer length (k)
 * k is limited to 2^8-1
 * @param numThreads
 */
void MantisSski::buildMPHF(uint32_t numThreads) {
    std::cerr << " storage size: " << contigSeq.size() << "\n";
    std::cerr << " contig size: " << contigStartIdx.size() << "\n";
    ContigKmerIterator kb(&contigSeq, &contigStartIdx, static_cast<uint8_t>(k), 0);
    ContigKmerIterator ke(&contigSeq, &contigStartIdx, static_cast<uint8_t>(k), contigSeq.size() - k + 1);
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
    ContigKmerIterator kb(&contigSeq, &contigStartIdx, static_cast<uint8_t>(k), 0);
    ContigKmerIterator ke(&contigSeq, &contigStartIdx, static_cast<uint8_t>(k), contigSeq.size() - k + 1);

    uint64_t contigIdx{0};
    // For every valid k-mer (i.e. every contig)
    uint32_t cntr{0}, prevPrcnt{0}, curPrcnt{0};
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
            std::cerr << "\r" << curPrcnt << "%";
            prevPrcnt = curPrcnt;
        }
    }

    // store prefixArr
    sdsl::store_to_file(prefixArr, outdir + "/prefix.bin");
    console->info("stored prefixArray.");
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