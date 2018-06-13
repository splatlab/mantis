#ifndef __COMPRESSED_SET_BIT_H__
#define __COMPRESSED_SET_BIT_H__
#include "codecfactory.h"
#include "intersection.h"

#include <vector>
#include <algorithm> /* for sort, random_shuffle */

using namespace SIMDCompressionLib;

template <typename IndexSizeT>
class CompressedSetBit {
public:
  CompressedSetBit(std::vector<uint32_t> idxList) {
    nSetInts = idxList.size();

    std::sort(idxList.begin(), idxList.end());
    //std::cout << "\nAfter sorting:\n";
    //for (size_t i = 0; i < idxList.size(); i++) std::cout << idxList[i] << " ";
    //std::cout << "\n";
    // We pick a CODEC
    //IntegerCODEC &codec = *CODECFactory::getFromName("s4-bp128-dm");
    IntegerCODEC &codec = *CODECFactory::getFromName("s4-fastpfor-d1");

    std::vector<uint32_t> dat(idxList.size() + 1024);
    //deltaCompressedSetList.resize(idxList.size() + 1024);
    size_t compressedsize = dat.size();
    codec.encodeArray(idxList.data(), idxList.size(), dat.data(),
                      compressedsize);

    dat.resize(compressedsize);
    dat.shrink_to_fit();
    nCompressedInts = compressedsize;
    data_.reset(new uint32_t[compressedsize]);
    std::copy(dat.begin(), dat.end(), data_.get());
  }

  void uncompress(std::vector<uint32_t>& idxList) {
    //IntegerCODEC &codec = *CODECFactory::getFromName("s4-bp128-dm");
    IntegerCODEC &codec = *CODECFactory::getFromName("s4-fastpfor-d1");
    idxList.resize(nSetInts);
    size_t compressedSize = nCompressedInts;//deltaCompressedSetList.size();
    size_t originalSize = nSetInts;
    codec.decodeArray(&data_[0], compressedSize,
                      idxList.data(), originalSize);
    idxList.resize(originalSize); //?? why do we need this?
  }

  size_t size_in_bytes() {
    //std::cout << "  size: " << deltaCompressedSetList.size() << " ";
    return sizeof(*this) + (sizeof(uint32_t) * nCompressedInts);
  }

  bool serialize(std::ostream& output) {
    output.write(reinterpret_cast<const char*>(&nSetInts), sizeof(nSetInts));
    output.write(reinterpret_cast<const char*>(&nCompressedInts), sizeof(nCompressedInts));
    output.write(reinterpret_cast<const char*>(&data_[0]), sizeof(uint32_t) * nCompressedInts);
    return true;
  }

  bool deserialize(std::istream& input) {
    input.read(reinterpret_cast<char*>(&nSetInts), sizeof(nSetInts));
    input.read(reinterpret_cast<char*>(&nCompressedInts), sizeof(nCompressedInts));
    data_.reset(new uint32_t[nCompressedInts]);
    input.read(reinterpret_cast<char*>(&data_[0]), sizeof(uint32_t) * nCompressedInts);
    return true;
  }

private:
  template <typename T>
  friend inline bool operator==(const CompressedSetBit<T>& lhs, const CompressedSetBit<T>& rhs);
  template <typename T>
  friend inline bool operator!=(const CompressedSetBit<T>& lhs, const CompressedSetBit<T>& rhs);
  //std::vector<uint32_t> deltaCompressedSetList;
  std::unique_ptr<uint32_t[]> data_;
  IndexSizeT nSetInts;
  IndexSizeT nCompressedInts;
};

template <typename T>
inline bool operator==(const CompressedSetBit<T>& lhs, const CompressedSetBit<T>& rhs) {
  return (lhs.nCompressedInts == rhs.nCompressedInts) ? (std::memcmp(&lhs.data_[0], &rhs.data_[0], lhs.nCompressedInts) == 0) : false;
}

template <typename T>
inline bool operator!=(const CompressedSetBit<T>& lhs, const CompressedSetBit<T>& rhs) {
  return !(lhs == rhs);
}

#endif
