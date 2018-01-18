#include <vector>
#include <cmath>
#include <unordered_map>
#include <iostream>
#include <fstream>
#include "sparsepp/spp.h"
#include "hashutil.h"
#include "bitvector.h"
int main(int argc, char *argv[])
{

  uint16_t num_samples = 2586;
  std::vector<uint64_t> queryCnt;
  std::string filename = argv[1];
  std::string outfileAddr = argv[2];
  size_t bootstrapCnt = std::stoull(argv[3]);
  BitVectorRRR eqcls(filename);
  size_t totalEqClsCnt = eqcls.bit_size() / num_samples; //222584822;
  std::cerr << "Total bit size: " << eqcls.bit_size() << "\n"
						<< "Total # of equivalence classes: " << totalEqClsCnt << "\n";
 
  size_t randCnt = 1000, blockSize = 100; // set it so that these two are divisible
  size_t numOfBlocks = randCnt/blockSize;
  size_t eqFreq[bootstrapCnt][numOfBlocks];
  sdslhash<sdsl::bit_vector> hasher;
  
  
  // bootstrap loops start here
  for (size_t bootstrapCntr = 0; bootstrapCntr < bootstrapCnt; bootstrapCntr++) {
		// reserve the space for each of the vectors for different block size hashes
		std::vector<std::vector<uint64_t>> blocks(numOfBlocks);
		for (int i = 0; i < numOfBlocks; i++) {
			blocks[i].reserve(totalEqClsCnt);
		}
		// choose randCnt(1000) random numbers between 0 and num_samples
		std::unordered_map<size_t, size_t> randIdx;
		std::vector<size_t> orderedRandIdx(randCnt);
		size_t i = 0;
		while (i < randCnt) {
			size_t currIdx = rand() % num_samples;
			if (randIdx.find(currIdx) != randIdx.end()) continue;
			randIdx[currIdx] = i;
			orderedRandIdx[i] = currIdx;
			i++;
		}
		std::sort(orderedRandIdx.begin(), orderedRandIdx.end());
		for (auto i : orderedRandIdx) {
			std::cerr << i << " ";
		}
		std::cerr << "\n";
		

		for (size_t eqclsCntr = 0; eqclsCntr < totalEqClsCnt; eqclsCntr++)
		{
			std::vector<sdsl::bit_vector> bvs(numOfBlocks);
			
			for (int i=0, siz = blockSize; siz <= randCnt; i++, siz+=blockSize) {
				sdsl::util::assign(bvs[i], sdsl::bit_vector(siz, 0));
			}
			//std::cerr << "\n" << eqclsCntr << " ";
			size_t i = 0;
			std::vector<size_t>::iterator idxPtr = orderedRandIdx.begin();
			while (i < num_samples && idxPtr != orderedRandIdx.end())
			{
				size_t bitCnt = std::min(num_samples, (uint16_t)64);
				//do get_int only if you have an index in that range, otherwise just skip and increase i
				if (idxPtr != orderedRandIdx.end() && (*idxPtr) < i+bitCnt) {
					//std::cerr << i << ":" ;
					size_t wrd = eqcls.get_int(eqclsCntr * num_samples + i, bitCnt);
					// fetch the bit value for all indices in that 64 bits of the fetched word
					while (idxPtr != orderedRandIdx.end() && (*idxPtr) < i+bitCnt) {	
						//std::cerr << *idxPtr << " -- ";
						uint8_t isOne = wrd >> ((*idxPtr)-i) & 0x01;
						if (isOne) {
							size_t idx = randIdx[(*idxPtr)];
							for (sdsl::bit_vector& bv : bvs) {
								if (idx < bv.size()) {
									//std::cerr << " s" << bv.size() << "," << idx << " - ";
									bv[idx] = 1;
								}
							}
						}
						idxPtr++;
					}
				}
				i += bitCnt;
			}
			
			/* for (int j = 0; j < bvs[4].size(); j++) {
				std::cerr << bvs[4][j]?'1':'0';
			}
			std::cerr << " h:" << hasher.getHash(bvs[4], bvs[4].capacity()/8);
 */			
		  for (int cntr = 0; cntr < blocks.size(); cntr++) {
				blocks[cntr].push_back(hasher.getHash(bvs[cntr], bvs[cntr].capacity()/8));
			}
			if (eqclsCntr != 0 && eqclsCntr % 1000000 == 0)
			{
				std::cout << "Total number of eqClses processed : " << eqclsCntr << "\n";
			}
		}
		//std::cerr << "\n\nblock size: " << blocks[0].size() << "\n";
		for (int cntr = 0; cntr < numOfBlocks; cntr++) {
			/* std::cerr << cntr <<" :\nbefore sort\n";
			for (auto v : blocks[cntr]) {
				std::cerr << v << " ";
			} */
			std::sort(blocks[cntr].begin(), blocks[cntr].end());
			/* std::cerr  << "\n" << cntr << " : after sort\n";
			for (auto v : blocks[cntr]) {
				std::cerr << v << " ";
			} */
			auto last = std::unique(blocks[cntr].begin(), blocks[cntr].end());
			/* std::cerr << "\n" << cntr << " : after unique size is : " << last - blocks[cntr].begin() << "\n";
			for (auto it = blocks[cntr].begin(); it != last; it++) {
				std::cerr << (*it) << " ";
			}
			std::cerr << "\n"; */
			eqFreq[bootstrapCntr][cntr] = last - blocks[cntr].begin();
		}
  }		  
  // end of the bootstrap loop

  std::ofstream outfile(outfileAddr);
  for (int bsTrialCntr = 0; bsTrialCntr < bootstrapCnt; bsTrialCntr++) {
  	for (int point = 0; point < numOfBlocks; point++) {
			outfile << eqFreq[bsTrialCntr][point] << "\t";
		}
		outfile << "\n";
  }
  outfile.close();
}
