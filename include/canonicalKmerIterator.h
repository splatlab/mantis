//
// Created by Fatemeh Almodaresi on 2019-03-14.
//

#ifndef MANTIS_CANONICALKMERITERATOR_H
#define MANTIS_CANONICALKMERITERATOR_H

#include "combine_canonicalKmer.hpp"
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

    ContigKmerIterator(sdsl::int_vector<2>* storage, sdsl::int_vector<>* offset,
                       uint8_t k, uint64_t startAt)
            : storage_(storage), offset_(offset), k_(k), curr_(startAt) {
//        std::cerr << "In constructor: " << curr_ << "\n";
        if (curr_ + k_ <= storage_->size()) {
            //nextValidPosition_();
            mer_.fromNum(storage_->get_int(2 * curr_, 2 * k_));
            if (contigCntr+1 < (*offset_).size()) {
                nextOffset_ = (*offset_)[contigCntr+1];
            }
        }
        // rcMer_ = mer_.get_reverse_complement();
    }

    bool advanceToValid() {
        nextValidPosition_();
        return (curr_ + k_ <= storage_->size());
    }

    ContigKmerIterator&
    operator=(ContigKmerIterator& other) { //}= default;//(sdsl::int_vector<>*
        // storage, sdsl::bit_vector* rank,
        // uint8_t k, uint64_t startAt) :
//        std::cerr << "=\n";
        storage_ = other.storage_;
        offset_ = other.offset_;
        k_ = other.k_;
        curr_ = other.curr_;
        mer_ = other.mer_;
        // rcMer_ = other.rcMer_;
        word_ = other.word_;
        contigCntr = other.contigCntr;
        nextOffset_ = other.nextOffset_;
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
//        std::cerr << "* " << curr_ << "\n";
        // word_ = (mer_.word(0) < rcMer_.word(0)) ? mer_.word(0) : rcMer_.word(0);
        word_ = mer_.getCanonicalWord();
        //std::cerr << mer_.to_str() << "\n";
        return word_;
    }

    difference_type pos() { return curr_; }

    bool isCanonical(){
        return mer_.fwWord() == mer_.getCanonicalWord() ;
    }

    bool isEndKmer() {
        size_t endPos = curr_ + k_ - 1;
        return endPos == nextOffset_;
    }

    pointer operator->() {
        word_ = mer_.getCanonicalWord(); //(mer_.word(0) < rcMer_.word(0)) ?
        // mer_.word(0) : rcMer_.word(0);
        return &word_;
    }
    bool operator==(const self_type& rhs) { return curr_ == rhs.curr_; }

    bool operator!=(const self_type& rhs) { return curr_ != rhs.curr_; }

    bool operator<(const self_type& rhs) { return curr_ < rhs.curr_; }

    bool operator<=(const self_type& rhs) { return curr_ <= rhs.curr_; }

    inline bool isValid() {
        size_t endPos = curr_ + k_ - 1;
        return !(endPos + 1 >= storage_->size() or endPos == nextOffset_);
    }

    inline bool isInvalid() {
        return !isValid();
    }

private:

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

    void advance_() {
        size_t endPos = curr_ + k_ - 1;
        //std::cerr << "advance_ " << endPos << " " << nextOffset_ << "\n";
        // If end of current valid kmer is the end of the contig (nextOffset_-1),
        // then next valid kmer is the first kmer in next contig
        if (endPos + 1 < storage_->size() and endPos == nextOffset_-1) {
            curr_ += k_;
            mer_.fromNum(storage_->get_int(2 * curr_, 2 * k_));
            contigCntr++;
            nextOffset_ = contigCntr+1 == offset_->size()?storage_->size():(*offset_)[contigCntr+1];
        } else {
            if (curr_ + k_ < storage_->size()) {
                int c = (*storage_)[curr_ + k_];
                mer_.shiftFw(c);
                ++curr_;
            } else {
                mer_.fromNum(storage_->get_int(2 * (storage_->size() - k_), 2 * k_));
                //std::cerr << "3before: " << curr_ << " ";
                curr_ = storage_->size() - k_ + 1;
               // std::cerr << "after: " << curr_ << "\n";

            }
        }
    }
    sdsl::int_vector<2>* storage_{nullptr};
    sdsl::int_vector<>* offset_{nullptr};
    uint8_t k_{0};
    uint64_t curr_{0};
    CanonicalKmer mer_;
    uint64_t word_{0};
    uint64_t contigCntr{0};
    uint64_t nextOffset_{0};
};

#endif //MANTIS_CANONICALKMERITERATOR_H
