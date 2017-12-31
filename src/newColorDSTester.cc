
//#include "compressedSetBit.h"
#include <vector>
#include <cmath>
#include <iostream>
#include <fstream>
#include "sparsepp/spp.h"
#include "hashutil.h"
#include "bitvector.h"

bool decodeAndValidate(std::map<size_t, BitVector*> &id2subPattern, std::vector<size_t> &blockList, 
BitVector& oldEqCls, size_t blockSize, size_t num_samples) {
  //decode
  BitVector bv(num_samples);
  size_t bvPos = 0;
  for (auto b : blockList) {
    //std::cout << b << ": ";
    BitVector* subbv = id2subPattern[b];
    size_t idx = 0;
    while (idx < blockSize) {  
      size_t bitCnt = std::min(blockSize - idx, (size_t)64);
      size_t wrd = subbv->get_int(idx, bitCnt);
      //std::cout << wrd << " " << bvPos << " " << std::min(num_samples-bvPos, bitCnt) << "\n";
      bv.set_int(bvPos, wrd, std::min(num_samples-bvPos, bitCnt));  
      idx += bitCnt;  
      bvPos += bitCnt;
    }
  }
  //compare
  size_t i = 0;
  size_t cntr = 0;
  while (i < num_samples) {
    size_t bitCnt = std::min(num_samples-i, (size_t)64);
    size_t wrd1 = bv.get_int(i, bitCnt);
    size_t wrd2 = oldEqCls.get_int(i, bitCnt);
    if (wrd1 != 0)
      std::cout << cntr << ": " << wrd1 << " " << wrd2 << "\n";
    if (wrd1 != wrd2) return false;
    i+= bitCnt;
    cntr++;
  }
  return true;
}
uint64_t addBlock(spp::sparse_hash_map<BitVector, std::pair<uint64_t, uint64_t>, sdslhash<BitVector>> &eqSubPattern_map, BitVector &subPattern)
{
  uint64_t eq_id = 0;
  //subPattern.calcHash();
  auto it = eqSubPattern_map.find(subPattern);
  // Find if the eqclass of the kmer is already there.
  // If it is there then increment the abundance.
  // Else create a new eq class.
  if (it == eqSubPattern_map.end()) { 
    eq_id = eqSubPattern_map.size();
    eqSubPattern_map.emplace(std::piecewise_construct, std::forward_as_tuple(subPattern), std::forward_as_tuple(eq_id, 1));
  }
  else { // eq class is seen before so increment the abundance.
    eq_id = it->second.first;
    it->second.second += 1; // update the abundance.
  }
  return eq_id;
}

/**
 * Calculates the disk/memory stats for flat split of color classes
 * two arguments
 * arg1: address to the eqClass file
 * arg2: subpattern (block) size
 * For example in mantis dir after make newColorDSTester run the following command:
 * ./newColorDSTester /mnt/scratch2/cdbg_cqf/raw_exact/eqclass_rrr.cls 256
 *
 **/
int main(int argc, char *argv[])
{

  uint16_t num_samples = 2586;
  std::vector<uint64_t> queryCnt;
  size_t totalQueryCnt = 0;
  std::string filename = argv[1];
  size_t subPatternLen = std::stoull(argv[2]);
  std::cerr << "SubPattern length: " << subPatternLen << "\n";
  BitVectorRRR eqcls(filename);
  size_t totalEqClsCnt = eqcls.bit_size() / num_samples; //222584822;
  std::cerr << "Total bit size: " << eqcls.bit_size() << "\ntotal # of equivalence classes: " << totalEqClsCnt << "\n";
  spp::sparse_hash_map<BitVector, std::pair<uint64_t, uint64_t>, sdslhash<BitVector>> subPatternMap;
  
  size_t validationCount = 1000000;
  for (size_t eqclsCntr = 0; eqclsCntr < totalEqClsCnt; eqclsCntr++)
  {
    size_t i = 0;
    std::map<uint64_t, uint8_t> idList;
    while (i < num_samples)
    {
      BitVector bv(subPatternLen);
      size_t bvPos = 0;
      size_t curSubPatternLen = std::min(num_samples - i, subPatternLen);
      size_t lastI = i + curSubPatternLen;
      while (i < lastI)
      {
        size_t bitCnt = std::min(lastI - i, (size_t)64);
        size_t wrd = eqcls.get_int(eqclsCntr * num_samples + i, bitCnt);
        bv.set_int(bvPos, wrd, bitCnt);
        /*for (size_t j = 0, bitCntr = bvPos; j < bitCnt; j++, bitCntr++) {
          if (wrd >> j & 0X01)
            bv.set(bitCntr);
            }*/
        bvPos += bitCnt;
        i += bitCnt;
      }
      uint64_t id = addBlock(subPatternMap, bv);
      if (idList.find(id) == idList.end())
      {
        idList[id] = 1;
      }
      else
        idList[id] += 1;
    }
    queryCnt.push_back(idList.size());
    totalQueryCnt += idList.size();
    if (eqclsCntr != 0 && eqclsCntr % 1000000 == 0)
    {
      std::cout << "\nTotal number of eqClses so far : " << eqclsCntr << "\nTotal average queries needed : " << totalQueryCnt / eqclsCntr << "\n";
    }
  }

  size_t subPatternBits = (ceil(log2(subPatternMap.size())));
  size_t subPatternCntPerEqCls = ceil(num_samples / (double)subPatternLen);
  size_t subPatternMapTotalSize = sizeof(subPatternMap) + subPatternMap.size() * (subPatternBits + 64); // key len + 2 size_t values that we store for each
  std::cerr << "Distribution of number of queries per eqCls:\n";
  std::map<uint64_t, uint64_t> queryDist;
  for (uint32_t qcnt : queryCnt)
  {
    if (queryDist.find(qcnt) != queryDist.end())
    {
      queryDist[qcnt]++;
    }
    else
      queryDist[qcnt] = 1;
  }
  for (auto it = queryDist.begin(); it != queryDist.end(); it++)
  {
    std::cerr << it->first << ": " << it->second << " , " << static_cast<int>(it->second * 100 / totalEqClsCnt) << "%\n";
  }
  size_t variableLenBits = 0;
  std::vector<std::pair<BitVector*, uint64_t>> subPatternFreqVec;
  subPatternFreqVec.reserve(subPatternMap.size());
  size_t sCntr = 0;
  
  // Insert subPatterns into a set sorted by their frequencies
  /* std::set<std::pair<decltype(subPatternMap::key_type), decltype(subPatternMap::value_type)>
  > subPatternSet_sortedByFreq(subPatternMap.begin(), subPatternMap.end(),
  [](std::pair<decltype(subPatternMap::key_type), decltype(subPatternMap::value_type)> a,
  std::pair<decltype(subPatternMap::key_type), decltype(subPatternMap::value_type)> b){
    return a.second > b.second;
  });
 */

  for (auto it = subPatternMap.begin(); it != subPatternMap.end(); it++) {
    subPatternFreqVec.emplace_back((BitVector*)(&(it->first)), it->second.second);
  }
  std::sort(subPatternFreqVec.begin(), subPatternFreqVec.end(),
            [](std::pair<const BitVector*, uint64_t> const a, std::pair<const BitVector*, uint64_t> const b) {
              return a.second > b.second;
            });
  
  // std::ofstream f;
  // f.open("subPatternFreq_" + std::to_string(subPatternLen) + ".dist");
  sCntr = 0;
  std::map<size_t, BitVector*> id2subPattern;
  for (auto it = subPatternFreqVec.begin(); it != subPatternFreqVec.end(); it++)
  {
    id2subPattern[sCntr] = it->first;
    subPatternMap[*(it->first)].first = sCntr; // update subPattern IDs
    variableLenBits += ceil(log2(sCntr + 1)) * it->second; 
    // f << sCntr++ << "," << *it << "\n";
    sCntr++;
  }
  // f.close();
  
  // validation section
  
/*   size_t blockCnt = ceil(num_samples/subPatternLen);
  for (size_t eqclsCntr = 0; eqclsCntr < validationCount; eqclsCntr++)
  {
    std::vector<uint64_t> newEqCls;
    newEqCls.reserve(blockCnt);
    size_t i = 0;
    BitVector oldEqCls(num_samples);
    size_t oldEqClsBitCntr = 0;
    size_t blckCntr = 0;
    while (i < num_samples)
    {
      BitVector bv(subPatternLen);
      size_t bvPos = 0;
      size_t curSubPatternLen = std::min(num_samples - i, subPatternLen);
      size_t lastI = i + curSubPatternLen;
      while (i < lastI)
      {
        size_t bitCnt = std::min(lastI - i, (size_t)64);
        size_t wrd = eqcls.get_int(eqclsCntr * num_samples + i, bitCnt);
        // std::cout << wrd << " ";
        bv.set_int(bvPos, wrd, bitCnt);
        for (size_t j = 0; j < bitCnt; j++) {
          if (wrd >> j & 0X01)
            oldEqCls.set(oldEqClsBitCntr);
          oldEqClsBitCntr++;
        }
        bvPos += bitCnt;
        i += bitCnt;
      }
      blckCntr++;
      newEqCls.push_back(subPatternMap[bv].first);
    }
    std::cout << "\n";
    if (!decodeAndValidate(id2subPattern, newEqCls, oldEqCls, subPatternLen, num_samples)) {
      std::cout << "Not Equal: " << eqclsCntr << "\n";
      std::exit(1);
    }
  }
  std::cout << "\nVALIDATION PASSED\n"; */
  size_t gigDenum = 8 * pow(1024, 3);
  size_t fixedLenTotal = subPatternMapTotalSize + subPatternBits * subPatternCntPerEqCls * totalEqClsCnt;
  size_t rsVarLenTotal = subPatternMapTotalSize + 2*variableLenBits;
  size_t huffmanVarLenTotal = subPatternMapTotalSize + 1.5*variableLenBits + subPatternMap.size()*64;
  std::cerr << "\n\n\nFinalResults:\n"
            << "SubPattern length: " << subPatternLen << "\nTotal number of experiments are: " << totalQueryCnt / totalEqClsCnt << "\nunique subPattern count: " << subPatternMap.size() << "\nsubpattern map size: " << subPatternMapTotalSize << "\n"
            << ceil(log2(subPatternMap.size())) << " bits instead of " << subPatternLen << " bits per each subPattern block."
            << "\nFixed length label: " << subPatternBits * subPatternCntPerEqCls << " bits per color class."
            << "\nVariable length label (on average): " << variableLenBits / totalEqClsCnt << " bits per color class."
            << "\n\nFixed length Overall: " << fixedLenTotal << " bits vs " << totalEqClsCnt * num_samples << " bits."
            << "or " << fixedLenTotal / gigDenum << "GB vs " << (totalEqClsCnt * num_samples) / gigDenum << "GB."
            << "\nSubpattern map total size: " << subPatternMapTotalSize << " bits or " << subPatternMapTotalSize/ gigDenum << "GB."
            << "\nRS Variable length Overall: " << rsVarLenTotal << " bits or " << (rsVarLenTotal) / gigDenum << "GB."
            << "\nHuffman Approximate Variable length Overall: " << huffmanVarLenTotal << " bits or " << (huffmanVarLenTotal) / gigDenum << "GB."
                                                                                                          "\n";
}
