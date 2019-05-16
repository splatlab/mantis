//
// Created by Fatemeh Almodaresi on 2019-03-14.
//

#ifndef MANTIS_CANONICALKMERITERATOR_H
#define MANTIS_CANONICALKMERITERATOR_H

#include "combine_canonicalKmer.hpp"

struct Sizes {
    uint32_t unitigSliceSize; // bps
    uint32_t bucketSize; // bits
    uint32_t slicePrefixSize; // bits
    void calcSlicePrefixSize() {
        slicePrefixSize =
                static_cast<uint32_t>(ceil(log2(static_cast<double>(unitigSliceSize))));
    }

    // don't change the 255 to 266.
    // 266 doesn't fit in 8 bits unless you consider 0 in 8-bit prefix as 1, 1 as 2, etc.
    Sizes(uint32_t bucketSizeIn = 1024, uint32_t unitigSliceSizeIn = 255) {
        bucketSize = bucketSizeIn * 8;
        unitigSliceSize = unitigSliceSizeIn;
        calcSlicePrefixSize();
        slicePrefixSize= static_cast<uint32_t>(ceil(log2(unitigSliceSize)));
        //unitigSliceSize+=slicePrefixSize;
    }
};

// adapted from :
// http://stackoverflow.com/questions/34875315/implementation-my-own-list-and-iterator-stl-c
class ContigKmerIterator {
public:
    using self_type = ContigKmerIterator;
    using value_type = uint64_t;
    using reference = value_type&;
    using pointer = value_type*;
    using iterator_category = std::forward_iterator_tag;
    using difference_type = int64_t;

    ContigKmerIterator(sdsl::int_vector<2>* storage, Sizes sizesIn,
                       uint8_t k, uint64_t startAt)
            : storage_(storage), sizes(sizesIn), k_(k), curr_(startAt) {
        CanonicalKmer::k(k_);
        if ((storage_->size() * 2) % sizes.bucketSize != 0) {
            std::cerr << "ERROR! Storage bits size is not divisible by bucketBits.\n";
            std::cerr << "Storage bits: " << storage_->size() * 2 <<
            " Bucket bits: " << sizes.bucketSize << "\n";
            std::exit(3);
        }
        //totalBuckets = storage_->size() * 2 / sizes.bucketSize;
//        std::cerr << "In constructor: " << curr_ << "\n";
        if (curr_ + k_ <= storage_->size()) {
            auto idxInBucket = curr_ % sizes.bucketSize;
            if (sizes.bucketSize - idxInBucket > sizes.slicePrefixSize) {
                sliceSize = storage_->get_int(2 * curr_, sizes.slicePrefixSize);
//                std::cerr << "constructor, curr: " << curr_ << " sliceSize: " << sliceSize << "\n";
                if (sliceSize) {
                    curr_ += (sizes.slicePrefixSize / 2); // move curr to the beginning of the contig
                    //nextValidPosition_();
                    mer_.fromNum(storage_->get_int(2 * curr_, 2 * k_));
                }
            }
        }
        // rcMer_ = mer_.get_reverse_complement();
    }

    /*bool advanceToValid() {
        nextValidPosition_();
        return (curr_ + k_ <= storage_->size());
    }
*/
    ContigKmerIterator&
    operator=(ContigKmerIterator& other) { //}= default;//(sdsl::int_vector<>*
        // storage, sdsl::bit_vector* rank,
        // uint8_t k, uint64_t startAt) :
//        std::cerr << "=\n";
        storage_ = other.storage_;
        k_ = other.k_;
        curr_ = other.curr_;
        mer_ = other.mer_;
        // rcMer_ = other.rcMer_;
        word_ = other.word_;
        //contigCntr = other.contigCntr;
        sliceSize = other.sliceSize;
        sizes = other.sizes;
        return *this;
    }

    ContigKmerIterator operator++() {
//        std::cerr << "++ " << curr_ << " ";
        ContigKmerIterator i = *this;
//        std::cerr << "i.curr_ " << i.curr_ << " ";
        advance_();
//        std::cerr << " -> " << curr_ << "\n";
        return i;
    }

    ContigKmerIterator operator++(int) {
//        std::cerr << "++int " << curr_ << " ";
        advance_();
//        std::cerr << curr_ << "\n";
        return *this;
    }

    reference operator*() {
        /*if (curr_-prev >= 100000000) {
            std::cerr << curr_ << " ";
            prev = curr_;
        }*/
//         word_ = (mer_.word(0) < rcMer_.word(0)) ? mer_.word(0) : rcMer_.word(0);
        word_ = mer_.getCanonicalWord();
//        std::cerr << curr_ << " seqsize: " << sliceSize << " " << mer_.to_str() << "\n";
        return word_;
    }

    difference_type pos() { return curr_; }

    bool isCanonical(){
        return mer_.fwWord() == mer_.getCanonicalWord() ;
    }

    /*bool isEndKmer() {
        size_t endPos = curr_ + k_ - 1;
        return endPos == nextOffset_;
    }
*/
    pointer operator->() {
        word_ = mer_.getCanonicalWord(); //(mer_.word(0) < rcMer_.word(0)) ?
        // mer_.word(0) : rcMer_.word(0);
        return &word_;
    }
    bool operator==(const self_type& rhs) { return curr_ == rhs.curr_; }

    bool operator!=(const self_type& rhs) { return curr_ != rhs.curr_; }

    bool operator<(const self_type& rhs) { return curr_ < rhs.curr_; }

    bool operator<=(const self_type& rhs) { return curr_ <= rhs.curr_; }

    /*inline bool isValid() {
        size_t endPos = curr_ + k_ - 1;
        return !(endPos + 1 >= storage_->size() or endPos == nextOffset_);
    }

    inline bool isInvalid() {
        return !isValid();
    }*/

    bool isSliceStart() { return !currInSlice_; }
    uint32_t getSliceSize() { return sliceSize; }
private:

/*
    void nextValidPosition_() {
        size_t endPos = curr_ + k_ - 1;
        if (endPos < storage_->size()) {
            // See if we cross a rank boundary
            bool crossesBoundary{false};
            bool isNextValid{false};
            size_t boundaryIndex{0};
            for (size_t i = curr_; i <= endPos; ++i) {
                if (i == nextOffset_) {
                    crossesBoundary = true;
                    boundaryIndex = i;
                    isNextValid = (i + k_) < storage_->size();
                    break;
                }
            }

            // If so, that's the start of the next valid k-mer
            // if that position is valid
            if (crossesBoundary) {
                if (isNextValid) {
                    //std::cerr << "1before: " << curr_ << " ";
                    curr_ = boundaryIndex; //+ 1; <- updated compared to rank
                    //std::cerr << "after: " << curr_ << "\n";
                    contigCntr++;
                    nextOffset_ = contigCntr+1 == offset_->size()?storage_->size():(*offset_)[contigCntr+1];
                } else {
                    // if that position is invalid, then go to the end.
                    goto endPos;
                }
            }
            // At this point, either curr_ points to the next valid
            // start position, or we have skipped over to the endPos label.
            mer_.fromNum(storage_->get_int(2 * curr_, 2 * k_));
            return;
        }

        endPos:
        // Fallthrough if we couldn't find a valid position.
        mer_.fromNum(storage_->get_int(2 * (storage_->size() - k_), 2 * k_));
        //std::cerr << "2before: " << curr_ << " ";
        curr_ = storage_->size() - k_ + 1;
        //std::cerr << "after: " << curr_ << "\n";

    }
*/
    void advance_() {
        // If we got to the end of the slice, jump to the next slice
        if (currInSlice_ + k_ == sliceSize) {
            curr_ += k_;
            auto idxInBucket = (curr_*2) % sizes.bucketSize;
            // if it is possible to have another slice in the bucket
            if (sizes.bucketSize - idxInBucket > sizes.slicePrefixSize) {
                sliceSize = storage_->get_int(2 * curr_, sizes.slicePrefixSize);
                if (!sliceSize) { // if there isn't a meaningful next slice (size = 0) go to next bucket
                    curr_ = (((curr_*2) / sizes.bucketSize) + 1) * sizes.bucketSize/2; // set curr to the start of next bucket
//                    std::cerr << "in if: " << curr_ << "\n";
                }

            } else { // go to next bucket
                curr_ = (((curr_*2) / sizes.bucketSize) + 1) * sizes.bucketSize/2; // set curr to the start of next bucket
//                std::cerr << "out if: " << curr_ << "\n";
            }
            if (curr_ >= storage_->size()) {
                curr_ = storage_->size() - k_ + 1;
                mer_.fromNum(storage_->get_int(2 * curr_, 2 * k_));
            } else { // if you're at a valid slice (either the same bucket or next bucket)
                sliceSize = storage_->get_int(2 * curr_, sizes.slicePrefixSize);
                curr_ += (sizes.slicePrefixSize / 2); // move curr to the beginning of the unitig
                currInSlice_ = 0;
                mer_.fromNum(storage_->get_int(2 * curr_, 2 * k_));
            }
        } else {
                int c = (*storage_)[curr_ + k_];
                mer_.shiftFw(c);
                ++curr_;
                ++currInSlice_;
        }
    }
    sdsl::int_vector<2>* storage_{nullptr};
    //sdsl::int_vector<>* offset_{nullptr};
    uint8_t k_{0};
    uint64_t curr_{0};
    uint64_t currInSlice_{0};
    CanonicalKmer mer_;
    uint64_t word_{0};
    //uint64_t bucketCnt{0};
    //uint64_t totalBuckets{0};
    //uint64_t nextOffset_{0};
    Sizes sizes;
    uint32_t sliceSize;
    uint64_t prev{0};
};

#endif //MANTIS_CANONICALKMERITERATOR_H
