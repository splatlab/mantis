
#include "compressedSetBit.h"
#include "bitvector.h"

int main(int argc, char* argv[]) {

  uint16_t num_samples = 2582;
  // set some values:

  std::set<uint32_t> randIdx;

  size_t cnt = 10;
  BitVector bv(num_samples);
  std::vector<uint32_t> idxList(cnt);
  std::cout << "original:\n";
  size_t i = 0;
  while (i < cnt) {
    uint32_t currIdx = rand() % num_samples;
    if (randIdx.find(currIdx) != randIdx.end()) continue;
    randIdx.insert(currIdx);
    bv.set(currIdx);
    idxList[i] = currIdx;
    std::cout << currIdx << " ";
    i++;
  }
  
  std::cout <<"\nidx size: " << idxList.size() << "\n";
  CompressedSetBit setBitList(idxList);

  vector<uint32_t> output;
  setBitList.uncompress(output);
  std::cout << "\nAfter compress&decompress size is " << output.size() << "\n";
  for (size_t i = 0; i < output.size(); i++) std::cout << output[i] << " ";
  std::cout << "\n";

}
