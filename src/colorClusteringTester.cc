
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
    eqSubPattern_map[subPattern] = std::make_pair(colorCls, 1);
    //eqSubPattern_map.emplace(std::piecewise_construct, std::forward_as_tuple(subPattern), std::forward_as_tuple(colorCls));
  }
  else { // eq class is seen before so increment the abundance.
    // second methodology : choose the eq. class with the maximum popularity
    //it->second.first = rand()%2?colorCls:it->second.first;
    it->second.second += 1; // update the abundance.
  }
  //return eq_id;
}

void addBlock(spp::sparse_hash_map<BitVector, std::vector<BitVector*>, sdslhash<BitVector>> &eqSubPattern_map, 
                  BitVector &subPattern,
                  BitVector* colorCls)
{
  auto it = eqSubPattern_map.find(subPattern);
  if (it == eqSubPattern_map.end()) { 
    std::vector<BitVector*> newVec;
    newVec.push_back(colorCls);
    eqSubPattern_map[subPattern] = newVec;
  }
  else { // eq class is seen before so increment the abundance.
    it->second.push_back(colorCls);
  }
}
// calculate hammingDistance between two bit vectors
size_t calcHammingDist(BitVector* b1, BitVector &b2, size_t num_samples, size_t minHammingDist = -1) {
  size_t hammingDistance = 0;
  size_t i = 0;
  while (i < num_samples) {
    size_t bitCnt = std::min(num_samples - i, (size_t)64);
    size_t wrd1 = (*b1).get_int(i, bitCnt);
    size_t wrd2 = b2.get_int(i, bitCnt);
    if (wrd1!=wrd2) {
      for (int j = 0; j < bitCnt; j++) {
        if ((wrd1 >> j & 0x01) != (wrd2 >> j & 0x01)) {
          hammingDistance++;
        }
      }
    }
    if (hammingDistance > minHammingDist) break;
    i += bitCnt;
  }
  return hammingDistance;
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
  std::cerr << "Top Most Frequent Eq. Clss Count: " << topFreqEqClsCnt << "\n";
  std::cerr << "loading file " << filename << "\n";
  BitVectorRRR eqcls(filename);
  size_t totalEqClsCnt = eqcls.bit_size() / num_samples; //222584822;
  std::cerr << "Total bit size: " << eqcls.bit_size() << "\ntotal # of equivalence classes: " << totalEqClsCnt << "\n";
  std::vector<spp::sparse_hash_map<BitVector, std::pair<BitVector*, uint64_t>, sdslhash<BitVector>>> subPatternMaps(ceil(subPatternCntPerEqCls));
  std::vector<spp::sparse_hash_map<BitVector, std::vector<BitVector*>, sdslhash<BitVector>>> subPatternMapsAndAllBvs(ceil(subPatternCntPerEqCls));
  //size_t validationCount = 1000000;
  // Go over the top most m frequent eq. classes. 
  // Record their IDs in a set for later use
  // Put each sub-pattern of them in its corresponding block map which points to the color eq. class bv
  std::ifstream eqClsFreq(eqclassFreq_filename);
  std::set<uint64_t> eqClsIDs;
  std::vector<BitVector> eqClss;
  eqClss.resize(topFreqEqClsCnt);
  for (size_t eqclsCntr = 0; eqclsCntr < topFreqEqClsCnt; eqclsCntr++)
  {
    std::string line;
    size_t eqClsID, tmp;
    eqClsFreq >> eqClsID >> tmp;
    eqClsID--; // in the file, the IDs start from 1 while eqCls indices are 0-based
    //std::cout << eqClsID << " ";
    eqClsIDs.insert(eqClsID);
    BitVector colorBv(num_samples);
    eqClss[eqclsCntr] = colorBv;
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
        eqClss[eqclsCntr].set_int(i, wrd, bitCnt);
        bvPos += bitCnt;
        i += bitCnt;
      }
      addBlock(subPatternMaps[blockCntr], bv, &eqClss[eqclsCntr]);
      addBlock(subPatternMapsAndAllBvs[blockCntr], bv, &eqClss[eqclsCntr]);
      blockCntr++;
    }
    // if (eqclsCntr != 0 && eqclsCntr % 1000000 == 0)
    //   std::cout << "\nTotal number of eqClses so far : " << eqclsCntr << "\n";
  }

   // calculate and store the distribution of collisions
  std::map<uint64_t, uint64_t> collisionMap;
  for (auto & subPatternMap : subPatternMaps) {
    for (auto it : subPatternMap) {
      if (collisionMap.find(it.second.second) == collisionMap.end()) {
        collisionMap[it.second.second] = 1;
      }
      else {
        collisionMap[it.second.second] += 1;
      }
    }
  }
  std::ofstream file("collisions-" + std::to_string(topFreqEqClsCnt) + ".dist");
  file << "CollisionCnt,BlockCnt\n";
  for (auto it : collisionMap) {
    file << it.first << "," << it.second << "\n";
  }
  file.close();
  //std::exit(1);

  // selecting 1000 random non-popular eq. classes
   size_t randCnt = 0; //1000;
/*  std::set<uint64_t> nonPopEqClsId;
  while (nonPopEqClsId.size() < randCnt) {
    size_t newID = rand() % totalEqClsCnt;
    if (eqClsIDs.find(newID) == eqClsIDs.end() && nonPopEqClsId.find(newID) == nonPopEqClsId.end()) {
      nonPopEqClsId.insert(newID);
    }
  }
 */
  std::cerr << "Looking for closest popular eq. class ...\n";
  // go over all eq classes other than the m most popular ones and find the popular vector with least distance
  size_t notFound = 0;
  std::map<size_t, size_t> hamDistMap;
  std::map<size_t, size_t> hamDist4AllMap;
  std::map<size_t, size_t> bruteForceHamDistMap;
  for (size_t eqclsCntr = 0; eqclsCntr < totalEqClsCnt; eqclsCntr++) {
  //for (size_t eqclsCntr : nonPopEqClsId) {
    //std::cout << eqclsCntr << "\n";
    if (eqClsIDs.find(eqclsCntr) == eqClsIDs.end()) {
      std::set<BitVector*> visitedBV;
      //std::cout << "not popular\n";      
      BitVector colorBv(num_samples);
      std::vector<BitVector> subPatterns;
      subPatterns.reserve(subPatternCntPerEqCls);
      size_t i = 0;
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
          colorBv.set_int(i, wrd, bitCnt);
          bvPos += bitCnt;
          i += bitCnt;
        }
        subPatterns.push_back(bv);
      }
      //std::cout << "Done with one color class\n";
      // find the block with minimum distance from the current color class and update stats for that distance
      size_t minHammingDist = -1;
      size_t minHammingDist4all = -1;
      bool found = false;
      for (size_t blockCntr=0; blockCntr<subPatterns.size(); blockCntr++) {
        auto it = subPatternMaps[blockCntr].find(subPatterns[blockCntr]);
        if (it != subPatternMaps[blockCntr].end()) {
          found = true;
          if (minHammingDist4all == 1) break;
          size_t hammingDist;
          // if (minHammingDist == 1 && minHammingDist4all == 1) break;

          // keeping the most popular
          /* hammingDist = calcHammingDist(it->second.first, colorBv, num_samples);
          if (hammingDist < minHammingDist) {
            minHammingDist = hammingDist;
          }
 */
          // keeping all the bitvectors:
          for (auto v : subPatternMapsAndAllBvs[blockCntr][subPatterns[blockCntr]]) {
            if (visitedBV.find(v) == visitedBV.end()) {
              visitedBV.insert(v);
              hammingDist = calcHammingDist(v, colorBv, num_samples, minHammingDist4all);
              if (hammingDist < minHammingDist4all) {
                minHammingDist4all = hammingDist;
              }
              if (minHammingDist4all == 1) break;
            }
          }     
        }
      }
      if (found) {
        // keeping the most popular
/*         if (hamDistMap.find(minHammingDist) != hamDistMap.end())
          hamDistMap[minHammingDist]++;
        else
          hamDistMap[minHammingDist] = 1;
 */
        // keeping all the bitvectors
        if (hamDist4AllMap.find(minHammingDist4all) != hamDist4AllMap.end())
          hamDist4AllMap[minHammingDist4all]++;
        else
          hamDist4AllMap[minHammingDist4all] = 1;
      }
      else { // if didn't find any similar blocks, increase notFound counter
          notFound++;
      }


      // brute-force search
/*       minHammingDist = -1;
      for (auto bv : eqClss) {
        size_t hammingDist = calcHammingDist(&bv, colorBv, num_samples);
        if (hammingDist < minHammingDist) {
          minHammingDist = hammingDist;
        }
        if (minHammingDist == 1) break;
      }
      if (bruteForceHamDistMap.find(minHammingDist) != bruteForceHamDistMap.end())
          bruteForceHamDistMap[minHammingDist]++;
        else
          bruteForceHamDistMap[minHammingDist] = 1; */
    }
    if (eqclsCntr != 0 && eqclsCntr % 1000000 == 0)  
      std::cerr << "\nTotal number of eqClses so far : " << eqclsCntr;
  }

  // write results for brute force
  std::ofstream bf_file("bf-"+std::to_string(topFreqEqClsCnt) + "-" + std::to_string(randCnt)+".dist");
  for (auto it : bruteForceHamDistMap) {
    bf_file << it.first << " " << it.second << "\n";
  }
  bf_file.close();

  // write results for keeping all bitvectors
  bf_file.open("keepAll-"+std::to_string(topFreqEqClsCnt) + "-" + std::to_string(randCnt)+".dist");
  for (auto it : hamDist4AllMap) {
    bf_file << it.first << " " << it.second << "\n";
  }
  bf_file.close();

// write results for keeping all bitvectors
  bf_file.open("keepMostPop-"+std::to_string(topFreqEqClsCnt) + "-" + std::to_string(randCnt)+".dist");
  for (auto it : hamDistMap) {
    bf_file << it.first << " " << it.second << "\n";
  }
  bf_file.close();

}