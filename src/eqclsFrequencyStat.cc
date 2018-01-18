#include <vector>
#include <cmath>
#include <unordered_map>
#include <iostream>
#include <fstream>
#include <time.h> 
#include "sparsepp/spp.h"
#include "hashutil.h"
#include "bitvector.h"
int main(int argc, char *argv[])
{

	srand(time(NULL));
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
  sdslhash<sdsl::bit_vector> hasher;
  
  size_t eqFreq[bootstrapCnt][numOfBlocks];
	std::vector<std::vector<std::vector<uint64_t>>> bsBlocks(bootstrapCnt);
	std::vector<std::unordered_map<uint16_t, uint16_t>> bsRandIdx(bootstrapCnt);
	std::vector<std::vector<uint16_t>> bsOrderedRandIdx(bootstrapCnt);
	std::vector<std::vector<sdsl::bit_vector>> bsbvs(bootstrapCnt);

  // bootstrap loops start here
  for (size_t bsCntr = 0; bsCntr < bootstrapCnt; bsCntr++) {
		// reserve the space for each of the vectors for different block size hashes
		bsBlocks[bsCntr].resize(numOfBlocks);
		for (int i = 0; i < numOfBlocks; i++) {
			bsBlocks[bsCntr][i].reserve(totalEqClsCnt);
		}
		// choose randCnt(1000) random numbers between 0 and num_samples
		bsOrderedRandIdx[bsCntr].resize(randCnt);
		size_t i = 0;
		while (i < randCnt) {
			size_t currIdx = rand() % num_samples;
			if (bsRandIdx[bsCntr].find(currIdx) != bsRandIdx[bsCntr].end()) continue;
			bsRandIdx[bsCntr][currIdx] = i;
			bsOrderedRandIdx[bsCntr][i] = currIdx;
			i++;
		}
		std::sort(bsOrderedRandIdx[bsCntr].begin(), bsOrderedRandIdx[bsCntr].end());
		for (auto i : bsOrderedRandIdx[bsCntr]) {
			std::cerr << i << " ";
		}
		std::cerr << "\n";
		bsbvs[bsCntr].resize(numOfBlocks);
		for (int i=0, siz = blockSize; siz <= randCnt; i++, siz+=blockSize) {
			sdsl::util::assign(bsbvs[bsCntr][i], sdsl::bit_vector(siz, 0));
		}
	}

	for (size_t eqclsCntr = 0; eqclsCntr < totalEqClsCnt; eqclsCntr++)
		{			
			//std::cerr << "\n" << eqclsCntr << " ";
			size_t i = 0;
			std::vector<std::vector<uint16_t>::iterator> bsIdxPtr(bootstrapCnt);
			for (int bsCntr=0; bsCntr<bootstrapCnt; bsCntr++) {
				bsIdxPtr[bsCntr] = bsOrderedRandIdx[bsCntr].begin();
			}
			while (i < num_samples)
			{
				size_t bitCnt = std::min(num_samples, (uint16_t)64);
				size_t wrd = eqcls.get_int(eqclsCntr * num_samples + i, bitCnt);
				// fetch the bit value for all indices in that 64 bits of the fetched word
				for (int bsCntr=0; bsCntr<bootstrapCnt; bsCntr++) {
					while (bsIdxPtr[bsCntr] != bsOrderedRandIdx[bsCntr].end() && (*bsIdxPtr[bsCntr]) < i+bitCnt) {	
						//std::cerr << *idxPtr << " -- ";
						uint8_t isOne = wrd >> ((*bsIdxPtr[bsCntr])-i) & 0x01;
						if (isOne) {
							size_t idx = bsRandIdx[bsCntr][(*bsIdxPtr[bsCntr])];
							for (sdsl::bit_vector& bv : bsbvs[bsCntr]) {
								if (idx < bv.size()) {
									//std::cerr << " s" << bv.size() << "," << idx << " - ";
									bv[idx] = 1;
								}
							}
						}
						bsIdxPtr[bsCntr]++;
					}
				}
				i += bitCnt;
			}
			
			/* for (int j = 0; j < bvs[4].size(); j++) {
				std::cerr << bvs[4][j]?'1':'0';
			}
			std::cerr << " h:" << hasher.getHash(bvs[4], bvs[4].capacity()/8);
 */		
			for (int bsCntr = 0; bsCntr < bootstrapCnt; bsCntr++) {	
				for (int cntr = 0; cntr < bsBlocks[bsCntr].size(); cntr++) {
					bsBlocks[bsCntr][cntr].push_back(hasher.getHash(bsbvs[bsCntr][cntr], bsbvs[bsCntr][cntr].capacity()/8));
				}
			}
			if (eqclsCntr != 0 && eqclsCntr % 1000000 == 0)
			{
				std::cout << "Total number of eqClses processed : " << eqclsCntr << "\n";
			}
		}
		for (int bsCntr = 0; bsCntr < bootstrapCnt; bsCntr++) {	
			
			//std::cerr << "\n\nblock size: " << blocks[0].size() << "\n";
			for (int cntr = 0; cntr < numOfBlocks; cntr++) {
				/* std::cerr << cntr <<" :\nbefore sort\n";
				for (auto v : blocks[cntr]) {
					std::cerr << v << " ";
				} */
				std::sort(bsBlocks[bsCntr][cntr].begin(), bsBlocks[bsCntr][cntr].end());
				/* std::cerr  << "\n" << cntr << " : after sort\n";
				for (auto v : blocks[cntr]) {
					std::cerr << v << " ";
				} */
				auto last = std::unique(bsBlocks[bsCntr][cntr].begin(), bsBlocks[bsCntr][cntr].end());
				/* std::cerr << "\n" << cntr << " : after unique size is : " << last - blocks[cntr].begin() << "\n";
				for (auto it = blocks[cntr].begin(); it != last; it++) {
					std::cerr << (*it) << " ";
				}
				std::cerr << "\n"; */
				eqFreq[bsCntr][cntr] = last - bsBlocks[bsCntr][cntr].begin();
			}
  	}		  
  // end of the bootstrap loop

	//output the results
  std::ofstream outfile(outfileAddr);
	for (int point = 0; point < numOfBlocks; point++) {
		outfile << bsbvs[0][point].size() << "\t";
	}
	outfile << "\n"; 
  for (int bsCntr = 0; bsCntr < bootstrapCnt; bsCntr++) {
  	for (int point = 0; point < numOfBlocks; point++) {
			outfile << eqFreq[bsCntr][point] << "\t";
		}
		outfile << "\n";
  }
  outfile.close();
}
