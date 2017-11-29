#ifndef __COMPRESSED_SET_BIT_H__
#define __COMPRESSED_SET_BIT_H__
#include "codecfactory.h"
#include "intersection.h"

#include <vector>
#include <algorithm> /* for sort, random_shuffle */

using namespace SIMDCompressionLib;

class CompressedSetBit {
public:
  CompressedSetBit(std::vector<uint32_t> idxList) {
    nSetBits = idxList.size();

    std::sort(idxList.begin(), idxList.end());
    //std::cout << "\nAfter sorting:\n";
    //for (size_t i = 0; i < idxList.size(); i++) std::cout << idxList[i] << " ";
    //std::cout << "\n";
    // We pick a CODEC
    IntegerCODEC &codec = *CODECFactory::getFromName("s4-bp128-dm");
    deltaCompressedSetList.resize(idxList.size() + 1024);
    size_t compressedsize = deltaCompressedSetList.size();
    codec.encodeArray(idxList.data(), idxList.size(), deltaCompressedSetList.data(),
                      compressedsize);
    deltaCompressedSetList.resize(compressedsize);
    deltaCompressedSetList.shrink_to_fit();
  }

  void uncompress(std::vector<uint32_t>& idxList) {
    IntegerCODEC &codec = *CODECFactory::getFromName("s4-bp128-dm");
    idxList.resize(nSetBits);
    size_t compressedSize = deltaCompressedSetList.size();
    size_t originalSize = nSetBits;
    codec.decodeArray(deltaCompressedSetList.data(), compressedSize,
                      idxList.data(), originalSize);
    idxList.resize(originalSize); //?? why do we need this?
  }

  size_t size_in_bytes() {
    //std::cout << "  size: " << deltaCompressedSetList.size() << " ";
    return sizeof(std::vector<uint32_t>) + (sizeof(uint32_t) * deltaCompressedSetList.size()); }
  bool serialize(std::ostream& output) {
    for (auto it = deltaCompressedSetList.begin(); it != deltaCompressedSetList.end(); it++)
      output << *it;
    output << "\n";
    return true;
  }

  bool deserialize(std::istream& input) {
    std::string val;
    input >> val;
    while (val != "\n") {
      uint32_t digVal(stoi(val));
      deltaCompressedSetList.push_back(digVal);
      input >> val;
    }
    return true;
  }
private:
  std::vector<uint32_t> deltaCompressedSetList;
  uint32_t nSetBits;
};

#endif
