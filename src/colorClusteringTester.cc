
//#include "compressedSetBit.h"
#include <vector>
#include <cmath>
#include <iostream>
#include <fstream>
#include <cstdlib>
#include "sparsepp/spp.h"
#include "hashutil.h"
#include "bitvector.h"


void addBlock(spp::sparse_hash_map<BitVector, std::pair<BitVector*, uint64_t>, sdslhash<BitVector>> &eqSubPattern_map, 
                  BitVector &subPattern,
                  BitVector* colorCls)
{
  auto it = eqSubPattern_map.find(subPattern);
  if (it == eqSubPattern_map.end()) { 
    eqSubPattern_map[subPattern] = std::make_pair(colorCls, eqSubPattern_map.size());
    //eqSubPattern_map.emplace(std::piecewise_construct, std::forward_as_tuple(subPattern), std::forward_as_tuple(colorCls));
  }
  else { // eq class is seen before so increment the abundance.
    it->second.first = rand()%2?colorCls:it->second.first;
    it->second.second += 1; // update the abundance.
  }
  //return eq_id;
}

double twoDecRound(double var) {
  size_t value = var * 100 + .5;
  return (double) value / 100.0;
}
/**
 * Calculates the disk/memory stats for hash-based clustering of color classes
 * two arguments
 * arg1: address to the directory containing color eq. classes and frequency file
 * arg2: subpattern (block) size
 * arg3 (optional, default=1M): count of top highest frequency eq. classes
 * For example in mantis dir after make colorClusteringTest run the following command:
 * ./colorClusteringTest /mnt/scratch2/cdbg_cqf/raw_exact 128 1000000
 *
 **/
int main(int argc, char *argv[])
{

  uint16_t num_samples = 2586;
  std::vector<uint64_t> queryCnt;
  std::string prefix = argv[1];
  std::string filename = prefix + "/eqclass_rrr.cls";
  std::string eqclassFreq_filename = prefix + "/tmpfile_sorted";
  size_t subPatternLen = std::stoull(argv[2]);
  size_t subPatternCntPerEqCls = ceil(num_samples / (double)subPatternLen);
  std::cerr << "SubPattern length: " << subPatternLen << "\n";
  size_t topFreqEqClsCnt = 1000000; //default
  if (argc > 3) {
    topFreqEqClsCnt = std::stoull(argv[3]);
  }
  std::cout<< "loading file " << filename << "\n";
  BitVectorRRR eqcls(filename);
  size_t totalEqClsCnt = eqcls.bit_size() / num_samples; //222584822;
  std::cerr << "Total bit size: " << eqcls.bit_size() << "\ntotal # of equivalence classes: " << totalEqClsCnt << "\n";
  std::vector<spp::sparse_hash_map<BitVector, std::pair<BitVector*, uint64_t>, sdslhash<BitVector>>> subPatternMaps(ceil(subPatternCntPerEqCls));
  //size_t validationCount = 1000000;
  // Go over the top most m frequent eq. classes. 
  // Record their IDs in a set for later use
  // Put each sub-pattern of them in its corresponding block map which points to the color eq. class bv
  std::ifstream eqClsFreq(eqclassFreq_filename);
  std::set<uint64_t> eqClsIDs;
  std::vector<BitVector> eqClss;
  eqClss.reserve(topFreqEqClsCnt);
  for (size_t eqclsCntr = 0; eqclsCntr < topFreqEqClsCnt; eqclsCntr++)
  {
    std::string line;
    size_t eqClsID, tmp;
    eqClsFreq >> eqClsID >> tmp;
    //std::cout << eqClsID << " ";
    eqClsIDs.insert(eqClsID);
    BitVector colorBv(num_samples);
    eqClss.push_back(colorBv);
    size_t i = 0;
    size_t blockCntr = 0;
    while (i < num_samples)
    {
      BitVector bv(subPatternLen);
      size_t bvPos = 0;
      size_t curSubPatternLen = std::min(num_samples - i, subPatternLen);
      size_t lastI = i + curSubPatternLen;
      while (i < lastI)
      {
        size_t bitCnt = std::min(lastI - i, (size_t)64);
        size_t wrd = eqcls.get_int(eqClsID * num_samples + i, bitCnt);
        bv.set_int(bvPos, wrd, bitCnt);
        eqClss.back().set_int(i, wrd, bitCnt);
        bvPos += bitCnt;
        i += bitCnt;
      }
      addBlock(subPatternMaps[blockCntr], bv, &eqClss.back());
      blockCntr++;
    }
    // if (eqclsCntr != 0 && eqclsCntr % 1000000 == 0)
    //   std::cout << "\nTotal number of eqClses so far : " << eqclsCntr << "\n";
  }
  std::cout << "Looking for closest popular eq. class ...\n";
  // go over all eq classes other than the m most popular ones and find the popular vector with least distance
  for (size_t eqclsCntr = 0; eqclsCntr < totalEqClsCnt; eqclsCntr++) {
    if (eqClsIDs.find(eqclsCntr) != eqClsIDs.end()) {
      
    }
  }
}