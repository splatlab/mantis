
#include "compressedSetBit.h"
#include "hashutil.h"
#include "bitvector.h"
#include "sdsl/bits.hpp"
#include "sparsepp/spp.h"

#include <time.h>
#include <random>
#include <sstream>
#include <algorithm>


void print_time_elapsed(string desc, struct timeval *start, struct timeval *end) {
    struct timeval elapsed;
    if (start->tv_usec > end->tv_usec) {
        end->tv_usec += 1000000;
        end->tv_sec--;
    }
    elapsed.tv_usec = end->tv_usec - start->tv_usec;
    elapsed.tv_sec = end->tv_sec - start->tv_sec;
    float time_elapsed = (elapsed.tv_sec * 1000000 + elapsed.tv_usec) / 1000000.f;
    std::cout << desc << "Total Time Elapsed: " << to_string(time_elapsed) << " seconds" << std::endl;
}

// @input
// cnt: setbit cnt 
// num_samples: total number of bits
// @return true if the output of delta_compression/decompression is the same as bitvector_rrr
//template<typename IndexSizeT>
bool validate(uint16_t cnt, uint16_t num_samples = 2586) {
    std::cout << "\n[Validate]\n";

    std::cout << "validate for " << cnt << " set bits out of " << num_samples << "\n";
    std::set<uint32_t> randIdx;

    BitVector bv(num_samples);
    std::vector<uint32_t> idxList(cnt);
    size_t i = 0;
    while (i < cnt) {
        //std::cout << "i:"<<i<<"\n";
        uint32_t currIdx = rand() % num_samples;
        //std::cout << "b" << currIdx << "\n";
        if (randIdx.find(currIdx) != randIdx.end()) continue;
        randIdx.insert(currIdx);
        bv.set(currIdx);
        idxList[i] = currIdx;
        //std::cout << currIdx << "\n";
        i++;
    }
    BitVectorRRR bvr(bv);
    std::vector<uint32_t> bvr_idxList;
    uint16_t wrdCnt = 64;
    for (uint16_t i = 0; i < num_samples; i += wrdCnt) {
        wrdCnt = std::min((uint16_t) 64, (uint16_t) (num_samples - i));
        uint64_t wrd = bvr.get_int(i, wrdCnt);
        for (uint16_t j = 0, idx = i; j < wrdCnt; j++, idx++) {
            if (wrd >> j & 0x01) {
                //std::cout << i << " " << j << "\n";
                bvr_idxList.push_back(idx);
            }
        }
    }
    //std::cout <<"\nidx size: " << idxList.size() << "\n";
    CompressedSetBit<uint32_t> setBitList(idxList);

    vector <uint32_t> output;
    setBitList.uncompress(output);
    //std::cout << "\nAfter compress&decompress size is " << output.size() << "\n";
    if (output.size() != bvr_idxList.size()) {
        std::cout << "rrr idx list size: " << bvr_idxList.size() << " deltac size: " << output.size() << "\n";
        return false;
    }
    for (size_t i = 0; i < output.size(); i++)
        if (output[i] != bvr_idxList[i]) {
            std::cout << i << " rrr idx: " << bvr_idxList[i] << " deltac: " << output[i] << "\n";
            return false;
        }
    return true;
}

/*
 * validates if the rrr bvs in the input directory all together have unique rows
 */
void validate_uniqueness(std::string &directory, size_t num_samples) {
    spp::sparse_hash_map<BitVector, uint64_t, sdslhash<BitVector>> eqMap;
    uint64_t cntr = 0;

    for (auto f = 0; f < 12; f++) {
        uint64_t curCntr = 0;
        std::string filename = directory + std::to_string(f) + "_eqclass_rrr.cls";
        BitVectorRRR eqcls(filename);
        //subPattern.calcHash();
        uint64_t bitcntr = 0;
        while (bitcntr < eqcls.bit_size()) {
            //std::cerr << "\n" << curCntr << "\n";

            BitVector eq(num_samples);
            size_t i = 0;
            while (i < num_samples) {
                size_t bitCnt = std::min(num_samples - i, (size_t) 64);
                size_t wrd = eqcls.get_int(curCntr * num_samples + i, bitCnt);
                for (size_t j = 0, curIdx = i; j < bitCnt; j++, curIdx++) {
                    if ((wrd >> j) & 0x01) {
                        eq.set(curIdx);

                    } /*else {
                        std::cerr << curIdx << " ";
                    }*/
                }
                i += bitCnt;
                bitcntr += bitCnt;
            }
            // Find if the eqclass of the kmer is already there.
            // If it is there then increment the abundance.
            // Else create a new eq class.
            if (eqMap.find(eq) == eqMap.end()) {
                eqMap.emplace(std::piecewise_construct, std::forward_as_tuple(eq), std::forward_as_tuple(cntr));
            } else { // eq class is seen before. Throw exception
                std::cerr << "\nFound an already existing equivalence class pattern. eq. number "
                          << cntr << " file " << f << " cur cntr: " << curCntr
                          << " before: " << eqMap[eq] << "\n";
                std::exit(1);
            }
            cntr++;
            curCntr++;
            if (curCntr % 1000000 == 0) {
                std::cerr << "passed " << curCntr << "th eq.\n";
            }
        }
    }
    std::cerr << "Validation passed!!! Found " << cntr << " unique eq classes.\n";


}

void compareCompressions(std::string &filename, size_t num_samples) {
    std::cout << "\n[CompareCompressions]\n";
    size_t gtBV = 0;
    size_t gtBVR = 0;
    size_t bvrGTbv = 0;
    size_t bvEQbvr = 0;
    size_t compressedSum = 0;
    size_t bvSum = 0;
    size_t bvrSum = 0;

    size_t roundIdxCnt = 0;
    size_t totalIdxCnt = 0;
    BitVectorRRR eqcls(filename);
    size_t totalEqClsCnt = eqcls.bit_size() / num_samples; //222584822;
    std::cout << "Total bit size: " << eqcls.bit_size() << "\ntotal # of equivalence classes: " << totalEqClsCnt
              << "\n";
    for (size_t eqclsCntr = 0; eqclsCntr < totalEqClsCnt; eqclsCntr++) {
        BitVector bv(num_samples);
        std::vector<uint32_t> idxList;
        idxList.reserve(num_samples);
        size_t i = 0;
        while (i < num_samples) {
            size_t bitCnt = std::min(num_samples - i, (size_t) 64);
            size_t wrd = eqcls.get_int(eqclsCntr * num_samples + i, bitCnt);
            for (size_t j = 0, curIdx = i; j < bitCnt; j++, curIdx++) {
                if ((wrd >> j) & 0x01) {
                    bv.set(curIdx);
                    idxList.push_back(curIdx);
                }
            }
            i += bitCnt;
        }
        //if (idxList.size() == 0) std::cerr << "Error!! " << eqclsCntr << " Shouldn't ever happen\n";
        totalIdxCnt += idxList.size();
        roundIdxCnt += idxList.size();
        if (eqclsCntr != 0 && (eqclsCntr) % 1000000 == 0) {
            std::cout << "\n\nTotal number of experiments are : " << eqclsCntr <<
                      "\nTotal average set bits: " << totalIdxCnt / eqclsCntr <<
                      "\nThis round average set bits: " << roundIdxCnt / 1000000 <<
                      "\nbv average size : " << bvSum / eqclsCntr <<
                      "\nbvr average size : " << bvrSum / eqclsCntr <<
                      "\ncompressed average size : " << compressedSum / eqclsCntr <<
                      "\ncompressed > bv : " << gtBV << " or " << (gtBV * 100) / eqclsCntr << "% of the time" <<
                      "\ncompressed > bvr : " << gtBVR << " or " << (gtBVR * 100) / eqclsCntr << "% of the time" <<
                      "\nbvr > bv : " << bvrGTbv << " or " << (bvrGTbv * 100) / eqclsCntr << "% of the time" <<
                      "\n";
            roundIdxCnt = 0;
        }

        BitVectorRRR bvr(bv);
        CompressedSetBit<uint64_t> setBitList(idxList);

        size_t bvSize = bv.bit_size()/8;//size_in_bytes();
        size_t bvrSize = bvr.bit_size()/8;//size_in_bytes();
        size_t compressedSize = setBitList.size_in_bytes();
        bvSum += bvSize;
        bvrSum += bvrSize;
        compressedSum += compressedSize;
        if (compressedSize > bvSize) gtBV++;
        if (compressedSize > bvrSize) gtBVR++;
        if (bvrSize > bvSize) bvrGTbv++;
        if (bvrSize == bvSize) bvEQbvr++;
        //vector<uint32_t> output;
        //setBitList.uncompress(output);
    }
    std::cout << "\n\n\nFinalResults:\n" <<
              "Total number of experiments are : " << totalEqClsCnt << "\n" <<
              "\nbv average size : " << bvSum / totalEqClsCnt <<
              "\nbvr average size : " << bvrSum / totalEqClsCnt <<
              "\ncompressed average size : " << compressedSum / totalEqClsCnt <<
              "\nHow many times compressed > bv : " << gtBV << " or " << (gtBV * 100) / totalEqClsCnt << "% of the time"
              <<
              "\nHow many times compressed > bvr : " << gtBVR << " or " << (gtBVR * 100) / totalEqClsCnt
              << "% of the time" <<
              "\nHow many times bvr > bv : " << bvrGTbv << " or " << (bvrGTbv * 100) / totalEqClsCnt << "% of the time"
              <<
              "\nHow many times bvr == bv : " << bvEQbvr << " or " << (bvEQbvr * 100) / totalEqClsCnt << "% of the time"
              <<
              "\n";
}

void compareCopies(size_t num_samples) {
    std::cout << "\n[CompareCopies]\n";

    struct timeval start, end;
    struct timezone tzp;

    size_t numOfCopies = 10;

    size_t bit_cnt = 1000000000;
    size_t i;
    sdsl::bit_vector bits(bit_cnt);
    srand(time(NULL));

    gettimeofday(&start, &tzp);
    for (auto i = 0; i < 1000000; i++) {
        bits[rand() % bit_cnt] = 1;
    }
    gettimeofday(&end, &tzp);
    print_time_elapsed("Set random bits done: ", &start, &end);

    sdsl::bit_vector bits2(bit_cnt);
    sdsl::bit_vector bits3(bit_cnt);
    sdsl::bit_vector bits4(bit_cnt + 1);

    gettimeofday(&start, &tzp);
    for (auto c = 0; c < numOfCopies; c++) {
        bits2 = bits;
    }
    gettimeofday(&end, &tzp);
    print_time_elapsed("MemCpy: ", &start, &end);

    gettimeofday(&start, &tzp);
    for (auto c = 0; c < numOfCopies; c++) {
        i = 0;
        while (i < bit_cnt) {
            if (bits[i])
                bits3[i] = 1;
            i++;
        }
    }

    gettimeofday(&end, &tzp);
    print_time_elapsed("Index by index!: ", &start, &end);

    gettimeofday(&start, &tzp);
    for (auto c = 0; c < numOfCopies; c++) {
        i = 0;
        size_t j = 1;
        while (i < bit_cnt) {
            size_t bitCnt = std::min(bit_cnt - i, (size_t) 64);
            size_t wrd = bits.get_int(i, bitCnt);
            bits4.set_int(j, wrd, bitCnt);
            i += bitCnt;
            j += bitCnt;
        }
    }
    gettimeofday(&end, &tzp);
    print_time_elapsed("Read int and write int: ", &start, &end);

}

void splitRows(std::string filename,
               std::string outdir,
               uint64_t num_samples = 2586) {
    //std::cerr << "1\n";
    std::string eqfile;
    std::ifstream eqlist(filename);
    std::ofstream bvlist(outdir+"/list.txt");
    uint64_t totalEqCls = 0;
    uint64_t rowcntr = 0;
    if (eqlist.is_open()) {
        while (getline(eqlist, eqfile)) {
            sdsl::rrr_vector<63> bvr;
            sdsl::load_from_file(bvr, eqfile);
            std::cerr << "loaded " << bvr.size() << "\n";
            totalEqCls = bvr.size()/num_samples;
            std::cerr << "total eq in this set: " << totalEqCls << "\t";
            std::cerr << "created\n";
            for (uint64_t eqcntr = 0; eqcntr < totalEqCls; eqcntr++) {
                sdsl::bit_vector eqcls(num_samples, 0);
                uint64_t i = 0;
                while (i < num_samples) {
                    uint64_t bitCnt = std::min(num_samples - i, (uint64_t) 64);
                    uint64_t wrd = bvr.get_int(eqcntr*num_samples+i, bitCnt);
                    //std::cerr << wrd << " ";
                    eqcls.set_int(i, wrd, bitCnt);
                    i+=bitCnt;
                }
                std::string outbvfile = outdir + "/row" + std::to_string(rowcntr) + ".sim.bf.bv";
                sdsl::store_to_file(eqcls, outbvfile);
                bvlist << outbvfile << "\n";
                rowcntr++;
            }
        }
        eqlist.close();
        bvlist.close();
    }

}

void splitColumns(std::string filename,
               std::string outdir,
               uint64_t totalEqClsCnt = 222000000,
               uint64_t num_samples = 2586) {
    std::string eqfile;
    std::ifstream eqlist(filename);
    std::ofstream bvlist(outdir+"/collist.txt");
    std::vector<sdsl::bit_vector> cols(num_samples);
    //std::vector<uint64_t> wrds(num_samples);
    for (auto &bv : cols) bv = sdsl::bit_vector(totalEqClsCnt, 0);
    if (eqlist.is_open()) {
        uint64_t accumTotalEqCls = 0, lineCntr{0};
        while (getline(eqlist, eqfile)) {
            //if (lineCntr >= 10) {
                sdsl::rrr_vector<63> bvr;
                std::cerr << "file: " << eqfile << "\n";
                sdsl::load_from_file(bvr, eqfile);
                std::cerr << "loaded " << bvr.size()
                          << " accumulative Eq Cnt: " << accumTotalEqCls << "\n";
                //std::cerr << "total eq in this set: " << totalEqCls << "\t";
                //std::cerr << "created\n";
                uint64_t b = 0;
                while (b < bvr.size() and accumTotalEqCls < totalEqClsCnt) {
                    uint64_t c{0}, bitcnt{0}, wrd{0}, sampleCntr{0};
                    while (c < num_samples) {
                        bitcnt = std::min(num_samples - c, (uint64_t) 64);
                        wrd = bvr.get_int(b + c, bitcnt);
                        for (uint64_t j = 0; j < bitcnt; j++) {
                            if ((wrd >> j) & 0x01) {
                                cols[sampleCntr][accumTotalEqCls] = 1;
                            }
                            sampleCntr++;
                        }
                        c += bitcnt;
                    }
                    /*for (uint64_t c = 0; c < num_samples; c++) {
                        if (bvr[(b + c)]) {
                            cols[c][accumTotalEqCls] = 1;
                        }
                    }*/
                    //std::cerr << "\n";
                    accumTotalEqCls++;
                    if (accumTotalEqCls % 1000000 == 0 || accumTotalEqCls > 222584820) {
                        std::cerr << accumTotalEqCls << " processed\n";
                    }
                    b += num_samples;
                }
            /*} else {
                lineCntr++;
                accumTotalEqCls += 20000000;
            }*/
        }
        std::cerr << "\n\nColumn bvs ready to be stored to file .. \n";
        uint64_t colcntr{0};
        for (auto& v : cols) {
            std::string outbvfile = outdir + "/col" + std::to_string(colcntr) + ".sim.bf.bv";
            sdsl::store_to_file(v, outbvfile);
            bvlist << outbvfile << "\n";
            colcntr++;
        }

        eqlist.close();
        bvlist.close();
    }

}


void decompress(std::string filename,
                std::string outfilename,
                uint64_t totalEqCls,
                uint64_t num_samples = 2586) {
    //std::cerr << "1\n";
    sdsl::rrr_vector<63> bvr;
    sdsl::load_from_file(bvr, filename);
    //sdsl::rank_support_rrr<1, 63> ranks(&bvr);
    //std::cerr << "1\n";
    //totalEqCls = 1;
    size_t prev = 0;
    size_t cur = 0;
    //std::cerr << "1\n";
    /*for (uint64_t i = 1; i < 21; i++) {
        cur = ranks(2586 * i);
        std::cout << "rank: " << cur << " " << cur - prev << "\n";
        prev = cur;
    }*/
    //std::exit(1);
    std::cerr << "loaded " << bvr.size() << "\n";
    sdsl::bit_vector eqcls(totalEqCls * num_samples, 0);
    std::cerr << "created\n";
    for (uint64_t i = 0; i < totalEqCls * num_samples; i += 64) {
        uint64_t bitCnt = std::min(totalEqCls * num_samples - i, (uint64_t) 64);
        uint64_t wrd = bvr.get_int(i, bitCnt);
        //std::cerr << wrd << " ";
        eqcls.set_int(i, wrd, bitCnt);
    }
    sdsl::store_to_file(eqcls, outfilename);
}

void compare_reordered_original(std::string filename1,
                                std::string filename2,
                                std::string order_filename,
                                uint64_t num_samples) {
    std::cout << "\n[compare_reordered_original]\n";

    struct timeval start, end;
    struct timezone tzp;

    uint64_t setBitCnt1{0}, setBitCnt2{0};
    std::ifstream orderFile(order_filename);
    std::vector<size_t> newOrder;
    while (!orderFile.eof()) {
        size_t val;
        orderFile >> val;
        //std::cout << val << " ";
        newOrder.push_back(val);
    }
    std::cout << "\n";
    sdsl::bit_vector eqcls;
    sdsl::bit_vector reordered_eqcls;
    sdsl::load_from_file(eqcls, filename1);
    sdsl::load_from_file(reordered_eqcls, filename2);
    size_t totalEqClsCnt = reordered_eqcls.size() / num_samples; //222584822;
    std::cout << "Total number of Eq. Clss: " << totalEqClsCnt << "\n";
    //totalEqClsCnt = 1000001;
    gettimeofday(&start, &tzp);

    for (size_t eqclsCntr = 0; eqclsCntr < totalEqClsCnt; eqclsCntr++) {

        std::vector<size_t> eqcls_setbits;
        std::vector<size_t> reordered_eqcls_setbits;
        eqcls_setbits.reserve(num_samples);
        reordered_eqcls_setbits.reserve(num_samples);

        // transpose the matrix of eq. classes
        size_t i = 0;
        while (i < num_samples) {
            size_t bitCnt = std::min(num_samples - i, (size_t) 64);
            size_t wrd1 = eqcls.get_int(eqclsCntr * num_samples + i, bitCnt);
            size_t wrd2 = reordered_eqcls.get_int(eqclsCntr * num_samples + i, bitCnt);
            for (size_t j = 0, curIdx = i; j < bitCnt; j++, curIdx++) {
                if ((wrd1 >> j) & 0x01) {
                    setBitCnt1++;
                    eqcls_setbits.push_back(newOrder[curIdx]);
                }
                if ((wrd2 >> j) & 0x01) {
                    setBitCnt2++;
                    reordered_eqcls_setbits.push_back(curIdx);
                }

            }
            i += bitCnt;
        }
        std::sort(eqcls_setbits.begin(), eqcls_setbits.end());
        std::sort(reordered_eqcls_setbits.begin(), reordered_eqcls_setbits.end());

        if (eqcls_setbits != reordered_eqcls_setbits) {
            std::cout << "ERROR! EQ CLSES in ROW: " << eqclsCntr << " NOT THE SAME.\n";
            std::exit(1);
        }
        if (eqclsCntr != 0 && (eqclsCntr) % 1000000 == 0) {
            gettimeofday(&end, &tzp);
            std::stringstream ss;
            ss << eqclsCntr << " eqclses processed, ";
            print_time_elapsed(ss.str(), &start, &end);
            gettimeofday(&start, &tzp);
        }

    }
    if (setBitCnt1 != setBitCnt2) {
        std::cout << "ERROR! EQ CLSES SET BITS NOT THE SAME: " << setBitCnt1 << " VS " << setBitCnt2 << "\n";
        std::exit(1);
    }
    std::cout << "\nValidation Passed!\n\n";
}

void reorder(std::string filename,
             std::string out_filename,
             std::vector<size_t> &newOrder,
             uint64_t num_samples) {
    std::cout << "\n[Reorder]\n";
    struct timeval start, end;
    struct timezone tzp;

    //std::cout << "\n";
    sdsl::bit_vector eqcls;
    sdsl::load_from_file(eqcls, filename);
    size_t totalEqClsCnt = eqcls.size() / num_samples; //222584822;
    std::cout << "Total number of Eq. Clss: " << totalEqClsCnt << "\n";
    //totalEqClsCnt = 1280000;
    sdsl::bit_vector reordered_eqcls(eqcls.size(), 0);
    gettimeofday(&start, &tzp);

    for (size_t eqclsCntr = 0; eqclsCntr < totalEqClsCnt; eqclsCntr++) {

        // transpose the matrix of eq. classes
        size_t i = 0;
        while (i < num_samples) {
            size_t bitCnt = std::min(num_samples - i, (size_t) 64);
            size_t wrd = eqcls.get_int(eqclsCntr * num_samples + i, bitCnt);
            for (size_t j = 0, curIdx = i; j < bitCnt; j++, curIdx++) {
                if ((wrd >> j) & 0x01) {
                    //std::cout << curIdx << ":" << newOrder[curIdx] << " ";
                    reordered_eqcls[eqclsCntr * num_samples + newOrder[curIdx]] = 1;
                }
            }
            i += bitCnt;
        }

        //std::cout << "\n";
        //std::exit(1);
        if (eqclsCntr != 0 && (eqclsCntr) % 1000000 == 0) {
            gettimeofday(&end, &tzp);
            std::stringstream ss;
            ss << eqclsCntr << " eqclses processed, ";
            print_time_elapsed(ss.str(), &start, &end);
            gettimeofday(&start, &tzp);
        }
    }
    std::cout << "\noutfilename: " << out_filename << "\n";
    for (auto idx = 0; idx < num_samples; idx++) {
        std::cout << reordered_eqcls[idx];
    }
    std::cout << "\n";
    for (size_t eqclsCntr = 0; eqclsCntr < 2; eqclsCntr++) {
        for (size_t i = 0; i < num_samples; i++) {
            std::cout << reordered_eqcls[num_samples * eqclsCntr + i];
        }
        std::cout << "\n";
    }
    sdsl::store_to_file(reordered_eqcls, out_filename);

}

void reorder(std::string filename,
             std::string out_filename,
             std::string order_filename,
             uint64_t num_samples) {

    std::ifstream orderFile(order_filename);
    std::vector<size_t> newOrder;
    newOrder.resize(num_samples);
    uint64_t cntr = 0;
    while (!orderFile.eof()) {
        size_t val;
        orderFile >> val;
        //std::cout << val << " ";
        newOrder[val] = cntr++;
        //newOrder.push_back(val);
    }

    reorder(filename, out_filename, newOrder, num_samples);
}

void build_fakeEq_And_randReord(std::string outfile,
                                uint64_t totalEqClsCnt = 1280000,
                                uint64_t num_samples = 2586,
                                uint64_t wordSize = 8) {
    std::cout << "\n[Build_fakeEq_And_randReord]\n";
    uint64_t bvSize = num_samples * totalEqClsCnt;
    sdsl::bit_vector eqcls(bvSize, 0);

    for (size_t eqclsCntr = 0; eqclsCntr < totalEqClsCnt; eqclsCntr++) {
        size_t i = 0;
        while (i < num_samples) {
            uint64_t wrd = rand() % 2 ? 0 : -1;
            size_t bitCnt = std::min(num_samples - i, (size_t) wordSize);
            eqcls.set_int(eqclsCntr * num_samples + i, wrd, bitCnt);
            //wrd = wrd?0:-1;
            i += bitCnt;
        }
    }

    for (size_t eqclsCntr = 0; eqclsCntr < 2; eqclsCntr++) {
        for (size_t i = 0; i < num_samples; i++) {
            std::cout << eqcls[num_samples * eqclsCntr + i];
        }
        std::cout << "\n";
    }
    sdsl::store_to_file(eqcls, outfile + ".eq");

    std::vector<uint64_t> randomizedIndex(num_samples);
    for (uint64_t i = 0; i < num_samples; i++)
        randomizedIndex[i] = i;
    std::random_shuffle(randomizedIndex.begin(), randomizedIndex.end());
    for (size_t i = 0; i < num_samples; i++) {
        std::cout << randomizedIndex[i] << " ";
    }
    std::cout << "\n";
    reorder(outfile + ".eq", outfile + "_randOrder.eq", randomizedIndex, num_samples);
}

void build_distMat(std::string filename,
                   std::string out_filename,
                   uint64_t num_samples) {
    std::cout << "\n[Build_distMat]\n";
    struct timeval start, end, row2colStart, colPairCompStart;
    struct timezone tzp;

    std::vector<std::vector<uint64_t>> hamDistMat(num_samples);
    for (auto &v : hamDistMat) {
        v.resize(num_samples);
        for (auto &d : v) d = 0;
    }

    std::vector<std::vector<double>> mutInfoMat(num_samples);
    std::vector<std::vector<double>> jointMat(num_samples);
    for (auto &v : mutInfoMat) {
        v.resize(num_samples);
        for (auto &d : v) d = 0;
    }
    for (auto &v : jointMat) {
        v.resize(num_samples);
        for (auto &d : v) d = 0;
    }

    //BitVectorRRR eqcls(filename);
    sdsl::bit_vector eqcls;
    sdsl::load_from_file(eqcls, filename);

    size_t totalEqClsCnt = eqcls.bit_size() / num_samples; //222584822;
    //totalEqClsCnt = 1280000;
    //sdsl::bit_vector uncompressed(totalEqClsCnt*num_samples, 0);
    std::cout << "Total bit size: " << eqcls.bit_size()
              << "\ntotal # of equivalence classes: " << totalEqClsCnt << "\n";

    std::vector<sdsl::bit_vector> cols(num_samples);
    for (auto &bv : cols) bv = sdsl::bit_vector(totalEqClsCnt, 0);

    gettimeofday(&start, &tzp);
    gettimeofday(&row2colStart, &tzp);

    for (size_t eqclsCntr = 0; eqclsCntr < totalEqClsCnt; eqclsCntr++) {

        // transpose the matrix of eq. classes
        size_t i = 0;
        while (i < num_samples) {
            size_t bitCnt = std::min(num_samples - i, (size_t) 64);
            size_t wrd = eqcls.get_int(eqclsCntr * num_samples + i, bitCnt);
            //uncompressed.set_int(eqclsCntr*num_samples+i, wrd, bitCnt);
            for (size_t j = 0, curIdx = i; j < bitCnt; j++, curIdx++) {
                if ((wrd >> j) & 0x01) {
                    cols[curIdx][eqclsCntr] = 1;
                }
            }
            i += bitCnt;
        }

        if (eqclsCntr != 0 && (eqclsCntr) % 1000000 == 0) {
            gettimeofday(&end, &tzp);
            std::stringstream ss;
            ss << eqclsCntr << " eqclses processed, ";
            print_time_elapsed(ss.str(), &start, &end);
            gettimeofday(&start, &tzp);
        }
    }
    //sdsl::store_to_file(uncompressed, "test.eq");
    gettimeofday(&end, &tzp);
    print_time_elapsed("DONE!, ", &row2colStart, &end);
    gettimeofday(&start, &tzp);

    gettimeofday(&colPairCompStart, &tzp);

    // calc ham distance
    size_t bitCnt = 64;
    for (size_t k = 0; k < num_samples; k++) {
        for (size_t j = k + 1; j < num_samples; j++) {
            uint64_t dist = 0;
            uint64_t mutprob[4];
            uint64_t denum = totalEqClsCnt;
            for (int b = 0; b < 4; b++) mutprob[b] = 0;
            uint64_t indivprob[2];
            for (int b = 0; b < 2; b++) indivprob[b] = 0;
            //calculate distance between column i and j
            size_t i = 0;
            while (i < totalEqClsCnt) {
                bitCnt = std::min(totalEqClsCnt - i, (size_t) 64);
                uint64_t wrd1 = cols[k].get_int(i, bitCnt);
                uint64_t wrd2 = cols[j].get_int(i, bitCnt);

                // compare words
                // xor & pop count
                dist += sdsl::bits::cnt(wrd1 ^ wrd2);

                indivprob[0] += sdsl::bits::cnt(wrd1);
                indivprob[1] += sdsl::bits::cnt(wrd2);
                mutprob[0] += sdsl::bits::cnt((~wrd1) & (~wrd2)) - (64 - bitCnt);
                mutprob[1] += sdsl::bits::cnt((~wrd1) & wrd2);
                mutprob[2] += sdsl::bits::cnt(wrd1 & (~wrd2));
                mutprob[3] += sdsl::bits::cnt(wrd1 & wrd2);
                i += bitCnt;
            }
            hamDistMat[k][j] = dist;
            hamDistMat[j][k] = dist;


            /* mutInfoMat[k][j] =
            (mutprob[0]?mutprob[0]*(log2(denum)+log2(mutprob[0])- log2(denum - indivprob[0])-log2(denum - indivprob[1])):0) +
            (mutprob[1]?mutprob[1]*(log2(denum)+log2(mutprob[1])- log2(denum - indivprob[0])-log2(indivprob[1])):0) +
            (mutprob[2]?mutprob[2]*(log2(denum)+log2(mutprob[2])- log2(indivprob[0])-log2(denum - indivprob[1])):0) +
            (mutprob[3]?mutprob[3]*(log2(denum)+log2(mutprob[3])- log2(indivprob[0])-log2(indivprob[1])):0);
            mutInfoMat[k][j] /= (double)(denum);
            mutInfoMat[j][k] = mutInfoMat[k][j]; */

            mutInfoMat[k][j] =
                    (mutprob[0] ? mutprob[0] * log2((double) (denum * mutprob[0]) /
                                                    (double) ((denum - indivprob[0]) * (denum - indivprob[1]))) : 0) +
                    (mutprob[1] ? mutprob[1] *
                                  log2((double) (denum * mutprob[1]) / (double) ((denum - indivprob[0]) * indivprob[1]))
                                : 0) +
                    (mutprob[2] ? mutprob[2] * log2((double) (denum * mutprob[2]) /
                                                    (double) ((indivprob[0]) * (denum - indivprob[1]))) : 0) +
                    (mutprob[3] ? mutprob[3] *
                                  log2((double) (denum * mutprob[3]) / (double) (indivprob[0] * (indivprob[1]))) : 0);
            mutInfoMat[k][j] /= (double) (denum);
            mutInfoMat[j][k] = mutInfoMat[k][j];

            jointMat[k][j] =
                    (mutprob[0] ? mutprob[0] * log2((double) (mutprob[0]) / (double) denum) : 0) +
                    (mutprob[1] ? mutprob[1] * log2((double) (mutprob[1]) / (double) denum) : 0) +
                    (mutprob[2] ? mutprob[2] * log2((double) (mutprob[2]) / (double) denum) : 0) +
                    (mutprob[3] ? mutprob[3] * log2((double) (mutprob[3]) / (double) denum) : 0);
            jointMat[k][j] /= (-1.0 * (double) (denum));
            jointMat[j][k] = jointMat[k][j];
            /*   std::cout << mutInfoMat[k][j] << " : "
                   << mutprob[0] << " "
                   << mutprob[1] << " "
                   << mutprob[2] << " "
                   << mutprob[3] << " "
                   << indivprob[0] << " "
                   << indivprob[1] << " "
                   << denum << "\n"; */

            //if (mutInfoMat[k][j]) std::cout << mutInfoMat[k][j] << "\n";
        }
        if (k % 100 == 0) {
            gettimeofday(&end, &tzp);
            std::stringstream ss;
            ss << k << ",";//<< j << ", ";
            print_time_elapsed(ss.str(), &start, &end);
            gettimeofday(&start, &tzp);
        }
    }
    gettimeofday(&end, &tzp);

    print_time_elapsed(" DONE!, ", &start, &end);
    gettimeofday(&colPairCompStart, &tzp);

    std::ofstream out(out_filename + ".ham");
    std::ofstream mutout(out_filename + ".mutinf"/* ".mutinfo" */);
    //std::ofstream jout(out_filename+".joint");
    //out << hamDistMat.size() << "\n";
    double mutInfoInv = 0;

    for (size_t i = 0; i < hamDistMat.size(); i++) {
        out << hamDistMat[i][0];
        mutInfoInv = (jointMat[i][0] - mutInfoMat[i][0]) < 0 ? 0 : (jointMat[i][0] - mutInfoMat[i][0]);
        mutout << mutInfoInv;
        for (size_t j = 1; j < hamDistMat[i].size(); j++) {
            out << "\t" << hamDistMat[i][j];
            mutInfoInv = (jointMat[i][j] - mutInfoMat[i][j]) < 0 ? 0 : (jointMat[i][j] - mutInfoMat[i][j]);
            mutout << "\t" << mutInfoInv;
            //jout << jointMat[i][j] << "\t";
            //if (hamDistMat[i][j] != 0) {

            //std::cout << j << " ";
            // out << i << "\t" << j << "\t" << hamDistMat[i][j] << "\n";
            //}
            //if (mutInfoMat[i][j] != 0) {
            // mutout <<  i << "\t" << j << "\t" << mutInfoMat[i][j] << "\t" << jointMat[i][j] << "\n";
            //}
        }
        out << "\n";
        mutout << "\n";
        //jout << "\n";
    }
}

void writeEq(std::string filename,
             std::string out_filename,
             uint64_t num_samples,
             uint64_t totalEqCls) {
    sdsl::rrr_vector<63> eqcls;
    sdsl::load_from_file(eqcls, filename);
    std::ofstream out(out_filename + ".matrix");
    size_t totalEqClsCnt = eqcls.size() / num_samples; //222584822;
    if (totalEqCls > 0) {
        totalEqClsCnt = totalEqCls;
    }
    for (size_t eqclsCntr = 0; eqclsCntr < totalEqClsCnt; eqclsCntr++) {
        size_t i = 0;
        while (i < num_samples) {
            size_t bitCnt = std::min(num_samples - i, (size_t) 64);
            size_t wrd = eqcls.get_int(eqclsCntr * num_samples + i, bitCnt);
            for (size_t j = 0; j < bitCnt; j++) {
                if (((wrd >> j) & 0x01))
                    out << "1";
                else
                    out <<" ";
                //out << ((wrd >> j) & 0x01) << "\t";
            }
            i += bitCnt;
        }
        out << "\n";
    }
    out.close();
}

int main(int argc, char *argv[]) {

    uint64_t num_samples = 2586;
    std::string command = argv[1];

    if (command == "validate") {
        std::cout << "validate for different number of set bits\n";
        for (uint16_t i = 1; i <= num_samples; i++) {
            if (!validate(i, num_samples)) {
                std::cerr << "ERROR: NOT VALIDATED\n";
                std::exit(1);
            }
        }
        std::cout
                << "SUCCESSFULLY VALIDATED THE COMPRESSION/DECOMPRESSION PROCESS.\n\n#####NEXT STEP#####\nCompare Sizes:\n";
    } else if (command == "decompress") {
        std::string filename = argv[2];
        std::string outfilename = argv[3];
        uint64_t totalEqCls = std::stoull(argv[4]);
        decompress(filename, outfilename, totalEqCls, num_samples);
    } else if (command == "compareCompressions") {
        std::string filename = argv[2];
        compareCompressions(filename, num_samples);
    } else if (command == "compareCopies")
        compareCopies(num_samples);
    else if (command == "distMat") {
        if (argc < 4) {
            std::cerr << "ERROR: MISSING LAST ARGUMENT\n";
            std::exit(1);
        }
        std::string filename = argv[2];
        std::string output_filename = argv[3];
        build_distMat(filename, output_filename, num_samples);
    } else if (command == "writeEq") {
        if (argc < 4) {
            std::cerr << "ERROR: MISSING LAST ARGUMENT\n";
            std::exit(1);
        }
        std::string filename = argv[2];
        std::string output_filename = argv[3];
        num_samples = stoull(argv[4]);
        uint64_t totalEqCls = 0;
        if (argc == 6)
            totalEqCls = std::stoull(argv[5]);

        writeEq(filename, output_filename, num_samples, totalEqCls);
    } else if (command == "reorder") {
        if (argc < 5) {
            std::cerr << "ERROR: MISSING ARGUMENT. 4 needed.\n";
            std::exit(1);
        }
        std::string filename = argv[2];
        std::string output_filename = argv[3];
        std::string order_filename = argv[4];
        reorder(filename, output_filename, order_filename, num_samples);
    } else if (command == "compare_reordered_original") {
        if (argc < 5) {
            std::cerr << "ERROR: MISSING ARGUMENT. 4 needed.\n";
            std::exit(1);
        }
        std::string filename1 = argv[2];
        std::string filename2 = argv[3];
        std::string order_filename = argv[4];
        compare_reordered_original(filename1, filename2, order_filename, num_samples);
    } else if (command == "fakeEq_randOrder") {
        std::string outfile = argv[2];
        uint64_t rows = 1048576;//2^20;
        build_fakeEq_And_randReord(outfile,
                                   rows,
                                   num_samples,
                                   8);
    } else if (command == "validate_uniqueness") {
        std::string dir = argv[2];
        num_samples = stoull(argv[3]);
        validate_uniqueness(dir, num_samples);
    } else if (command == "splitRows") {
        std::string filename = argv[2];
        std::string outdir = argv[3];
        num_samples = stoull(argv[4]);
        splitRows(filename, outdir, num_samples);
    } else if (command == "splitColumns") {
        std::string filename = argv[2];
        std::string outdir = argv[3];
        num_samples = stoull(argv[4]);
        uint64_t num_eqs = stoull(argv[5]);
        splitColumns(filename, outdir, num_eqs, num_samples);
    } else if (command == "quickvalidation") {
        std::string filename = argv[2];
        sdsl::bit_vector eqcls;
        sdsl::load_from_file(eqcls, filename);
        std::ofstream out(filename + ".seq");
        uint64_t i = 0;
        out << "start " << eqcls.bit_size() << "\n";
        while (i < eqcls.bit_size()) {
            size_t wrd = eqcls.get_int(i, 64);
            for (size_t j = 0; j < 64; j++) {
                if (((wrd >> j) & 0x01)) {
                    std::cerr << (i+j) << " ";
                }
            }
            i += 64;
        }
        out.close();

    } else {
        std::cerr << "ERROR: NO COMMANDS PROVIDED\n"
                  << "OPTIONS ARE: validate, compareCompressions, compareCopies, reorder\n";
        std::exit(1);
    }

}

