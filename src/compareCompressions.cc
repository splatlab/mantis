
#include "compressedSetBit.h"
#include "bitvector.h"

// @input
// cnt: setbit cnt 
// num_samples: total number of bits
// @return true if the output of delta_compression/decompression is the same as bitvector_rrr
//template<typename IndexSizeT>
bool validate(uint16_t cnt, uint16_t num_samples=2586) {

  std::cout << "validate for " << cnt << " set bits out of " << num_samples <<"\n";
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
  for (uint16_t i = 0; i < num_samples; i+=wrdCnt) {
	wrdCnt = std::min((uint16_t)64, (uint16_t)(num_samples - i));
    uint64_t wrd = bvr.get_int(i, wrdCnt);
    for (uint16_t j = 0, idx=i; j < wrdCnt; j++, idx++) {
        if (wrd >> j & 0x01) {
			//std::cout << i << " " << j << "\n";
				bvr_idxList.push_back(idx);
		}
    }  
  }
  //std::cout <<"\nidx size: " << idxList.size() << "\n";
  CompressedSetBit<uint32_t> setBitList(idxList);
  
  vector<uint32_t> output;
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


int main(int argc, char *argv[]) {


  uint16_t num_samples = 2586;
  std::cout << "validate for different number of set bits\n";
  for (uint16_t i = 1; i <= num_samples; i++){
  	 	if (!validate(i, num_samples)) {
				std::cerr << "ERROR: NOT VALIDATED\n";
				std::exit(1);
		}
  }
  std::cout << "SUCCESSFULLY VALIDATED THE COMPRESSION/DECOMPRESSION PROCESS.\n\n#####NEXT STEP#####\nCompare Sizes:\n";
  
  size_t gtBV = 0;
  size_t gtBVR = 0;
  size_t bvrGTbv = 0;
  size_t bvEQbvr = 0;
  size_t compressedSum = 0;
  size_t bvSum = 0;
  size_t bvrSum = 0;

  size_t roundIdxCnt = 0;
  size_t totalIdxCnt = 0;
  std::string filename = argv[1];
  BitVectorRRR eqcls(filename);
  size_t totalEqClsCnt = eqcls.bit_size()/num_samples; //222584822;
  std::cout << "Total bit size: " << eqcls.bit_size() << "\ntotal # of equivalence classes: " << totalEqClsCnt << "\n";
  for (size_t eqclsCntr = 0; eqclsCntr < totalEqClsCnt; eqclsCntr++) {
    BitVector bv(num_samples);
    std::vector<uint32_t> idxList;
    idxList.reserve(num_samples);
	size_t i = 0;
    while (i < num_samples) {
      size_t bitCnt = std::min(num_samples-i, (size_t)64);
      size_t wrd = eqcls.get_int(eqclsCntr*num_samples+i, bitCnt);
      for (size_t j = 0, curIdx = i; j < bitCnt; j++, curIdx++) {
        if ((wrd >> j) & 0x01) {
          bv.set(curIdx);
          idxList.push_back(curIdx);
        }
      }
      i+=bitCnt;
    }
    //if (idxList.size() == 0) std::cerr << "Error!! " << eqclsCntr << " Shouldn't ever happen\n";
    totalIdxCnt += idxList.size();
    roundIdxCnt += idxList.size();
    if (eqclsCntr != 0 && (eqclsCntr) % 1000000 == 0) {
      std::cout << "\n\nTotal number of experiments are : " << eqclsCntr <<
        "\nTotal average set bits: " << totalIdxCnt /eqclsCntr <<
        "\nThis round average set bits: " << roundIdxCnt/1000000 <<
                                                        "\nbv average size : " << bvSum/eqclsCntr <<
                                                        "\nbvr average size : " << bvrSum/eqclsCntr <<
                                                        "\ncompressed average size : " << compressedSum/eqclsCntr <<
                                                        "\ncompressed > bv : " << gtBV << " or " << (gtBV*100)/eqclsCntr << "% of the time" <<
                                                        "\ncompressed > bvr : " << gtBVR << " or " << (gtBVR*100)/eqclsCntr << "% of the time" <<
                                                        "\nbvr > bv : " << bvrGTbv << " or " << (bvrGTbv*100)/eqclsCntr << "% of the time" <<
                                                        "\n";
      roundIdxCnt = 0;
    }

    BitVectorRRR bvr(bv);
    CompressedSetBit<uint64_t> setBitList(idxList);

    size_t bvSize = bv.size_in_bytes();
    size_t bvrSize = bvr.size_in_bytes();
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
    "\nbv average size : " << bvSum/totalEqClsCnt <<
    "\nbvr average size : " << bvrSum/totalEqClsCnt <<
    "\ncompressed average size : " << compressedSum/totalEqClsCnt <<
    "\nHow many times compressed > bv : " << gtBV << " or " << (gtBV*100)/totalEqClsCnt << "% of the time" <<
    "\nHow many times compressed > bvr : " << gtBVR << " or " << (gtBVR*100)/totalEqClsCnt << "% of the time" <<
    "\nHow many times bvr > bv : " << bvrGTbv << " or " << (bvrGTbv*100)/totalEqClsCnt << "% of the time" <<
    "\nHow many times bvr == bv : " << bvEQbvr <<  " or " << (bvEQbvr*100)/totalEqClsCnt << "% of the time" <<
    "\n";

}
