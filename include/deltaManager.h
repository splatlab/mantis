//
// Created by Fatemeh Almodaresi on 8/17/18.
//

#ifndef MANTIS_DELTAMANAGER_H
#define MANTIS_DELTAMANAGER_H

#include <vector>
#include <math.h>
#include <iostream>
//#include <bitset>
#include <math.h>
#include <exception>
#include "sdsl/bit_vectors.hpp"

struct DeltaManagerException : public std::exception {
private:
    std::string message_;
public:

    DeltaManagerException(const std::string& message) : message_(message) {}
    const char * what () const throw () {
        return message_.c_str();
    }
};

class DeltaManager {
public:

    DeltaManager(uint64_t numSamplesIn,
                 uint64_t approximateColorClsCnt,
                 uint64_t approximateAvgDeltaCntPerColorCls) : numSamples(numSamplesIn) {
        slotWidth = ceil(log2(numSamples));
        if (slotWidth * (numSamples+1) < 64) {
            slotsPerColorCls = numSamples + 1;
        } else if (slotWidth * (approximateAvgDeltaCntPerColorCls+1) < 64) {
            slotsPerColorCls = (64 / slotWidth + 1) + 1;
        } else {
            slotsPerColorCls = approximateAvgDeltaCntPerColorCls + 1; // 1 for storing count of deltas per index
        }
        deltas.reserve(approximateColorClsCnt * slotsPerColorCls);
        // assumption: count of slots * their width is greater than 64 bits
        slotsPerColorClsWithPtrs = ((slotsPerColorCls - 1) * slotWidth - 64) / slotWidth + 1;
        colorCnt = 0;
    }

    void insertDeltas(uint64_t colorId, const std::vector<uint64_t> &dlta) {

        // see if we need to split deltas between main DS and heap
        // TODO not an expected behaviour, but we have no choice
        //if (colorId > colorCnt) throw DeltaManagerException("colorId > colorCnt");

        auto startBit = colorId*slotsPerColorCls*slotWidth;
        auto nextStartBit = (colorId+1)*slotsPerColorCls*slotWidth;
        // We always want to assign slotsPerColorCls slots to each index
        // even if deltas in that index are fewer than the avg num of deltas
        if (colorId >= colorCnt) {
            while (deltas.size() < nextStartBit/64+1) {
                deltas.push_back(0);
            }
            colorCnt = colorId+1;
        } else {
            // take care of deleting the pointer in case of previously creating one here:
            deletePtr(colorId);
        }

        // now assuming we've already reserved the space for the new colorId, insert deltas
        uint64_t mainDSDeltaCnt = dlta.size() < slotsPerColorCls ? dlta.size() : slotsPerColorClsWithPtrs - 1;
        insertValIntoDeltaV(startBit, dlta.size(), slotWidth); // insert count of deltas
        totDeltaCnt += dlta.size();
        startBit += slotWidth;
        for (uint64_t i = 0; i < mainDSDeltaCnt; i++) { // insert values into main DS
            insertValIntoDeltaV(startBit, dlta[i], slotWidth);
            startBit += slotWidth;
        }
        if (mainDSDeltaCnt < dlta.size()) { // in case count of deltas exceeds the reserved count
            // store the rest in an array in heap
            uint64_t *theRest = new uint64_t[(uint64_t)std::ceil((double)((dlta.size() - mainDSDeltaCnt) * slotWidth) / 64.0)]();
            insertValIntoDeltaV(startBit, reinterpret_cast<uint64_t>(theRest),
                                64); // store the pointer to the heap in the main DS
            insertValIntoHeap(dlta, mainDSDeltaCnt, theRest, slotWidth);
        }
    }

    std::vector<uint64_t> getDeltas(uint64_t colorId) {
        if (colorId > colorCnt) throw DeltaManagerException("colorId > colorCnt");
        std::vector<uint64_t> res;
        uint64_t startBit = colorId * slotWidth * slotsPerColorCls; // index for next color
        uint64_t deltaCnt = getValFromMDeltaV(startBit, slotWidth);
        uint64_t mainDSDeltaCnt = deltaCnt < slotsPerColorCls ? deltaCnt : slotsPerColorClsWithPtrs - 1;
        startBit += slotWidth;
        for (uint64_t i = 0; i < mainDSDeltaCnt; i++) { // insert values into main DS
            uint64_t delta = getValFromMDeltaV(startBit, slotWidth);
            res.push_back(delta);
            startBit += slotWidth;
        }
        if (mainDSDeltaCnt < deltaCnt) { // in case count of deltas exceeds the reserved count
            // fetch the pointer to the heap in the main DS
            uint64_t *theRestV = reinterpret_cast<uint64_t *>(getValFromMDeltaV(startBit, 64));
            getValFromHeap(res, deltaCnt - mainDSDeltaCnt, theRestV, slotWidth);
        }
        return res;
    }

    void swapDeltas(uint64_t colorId1, uint64_t colorId2) {
        std::vector<uint64_t> c1deltas = getDeltas(colorId1);
        std::vector<uint64_t> c2deltas = getDeltas(colorId2);
        insertDeltas(colorId1, c2deltas);
        insertDeltas(colorId2, c1deltas);
    }

    bool serialize(std::string prefix) {
        std::string deltabv_file = prefix + "/delta.bv";
        std::string boundarybv_file = prefix + "/boundary.bv";

        sdsl::int_vector<> deltabv(totDeltaCnt, 0, slotWidth);
        sdsl::bit_vector boundarybv(totDeltaCnt, 0);
        uint64_t j = 0;
        for (uint64_t i = 0; i < colorCnt; i++) {
            auto dltas = getDeltas(i);
            for (auto dlt : dltas) {
                deltabv[j] = dlt;
                j++;
            }
            boundarybv[j-1] = 1;
        }
        bool deltabvSuccessfullyStored = sdsl::store_to_file(deltabv, deltabv_file);
        bool boundarybvSuccessfullyStored = sdsl::store_to_file(boundarybv, boundarybv_file);

        return deltabvSuccessfullyStored and boundarybvSuccessfullyStored;
    }

private:
    uint64_t numSamples;
    uint64_t slotWidth;
    uint64_t slotsPerColorCls;
    uint64_t slotsPerColorClsWithPtrs;
    std::vector<uint64_t> deltas;
    uint64_t colorCnt;
    uint64_t totDeltaCnt{0};

    void deletePtr(uint64_t colorId) {
        uint64_t startBit = colorId * slotWidth * slotsPerColorCls; // index for next color
        uint64_t deltaCnt = getValFromMDeltaV(startBit, slotWidth); // get num of deltas in this slot
        totDeltaCnt -= deltaCnt;
        uint64_t mainDSDeltaCnt = deltaCnt < slotsPerColorCls ? deltaCnt : slotsPerColorClsWithPtrs - 1;
        if (mainDSDeltaCnt < deltaCnt) { // in case count of deltas exceeds the reserved count
            // fetch the pointer to the heap in the main DS and DELETE it
            startBit += slotWidth * (mainDSDeltaCnt + 1);
            uint64_t *theRestV = reinterpret_cast<uint64_t *>(getValFromMDeltaV(startBit, 64));
            delete theRestV;
        }
    }
    // width is limited to 64 (word size)
    bool insertValIntoDeltaV(uint64_t startBit, uint64_t val, uint64_t width) {
        uint64_t mask = width < 64? (((uint64_t)1 << width) - 1) : -1;
        if (startBit >= deltas.size() * 64) throw DeltaManagerException("startBit exceeds bit_size");
        if (width == slotWidth and val >= numSamples) {
            std::string msg = "delta index is larger than num_samples. val:"+
                    std::to_string(val)+
                              " num_samples:" +std::to_string(numSamples);
            throw DeltaManagerException(msg);
        }

        uint64_t &wrd = deltas[startBit / 64];
        uint64_t startBitInwrd = startBit % 64;
        wrd &= ~(mask << startBitInwrd);
        wrd |= (val << startBitInwrd);
        uint64_t shiftRight = 64 - startBitInwrd;
        if (shiftRight < width) {
            if (startBit >= (deltas.size() - 1) * 64) throw DeltaManagerException("startBit exceeds bit_size");
            auto &nextWrd = deltas[startBit/64+1];
            nextWrd &= ~(mask >> shiftRight);
            nextWrd |= (val >> shiftRight);
        }

        return true;
    }

    // width is limited to 64 (word size)
    uint64_t getValFromMDeltaV(uint64_t startBit, uint64_t width) {
        if (startBit >= deltas.size() * 64) throw DeltaManagerException("startBit exceeds bit_size");
        uint64_t mask = width < 64? (((uint64_t)1 << width) - 1) : -1;
        uint64_t res{0};
        uint64_t wrd = deltas[startBit / 64];
        uint64_t startBitInwrd = startBit % 64;
        res |= (mask & (wrd >> startBitInwrd));

        uint64_t shiftLeft = 64 - startBitInwrd;
        if (shiftLeft < width) {
            if (startBit >= (deltas.size() - 1) * 64) throw DeltaManagerException("startBit exceeds bit_size");
            wrd = deltas[startBit / 64 + 1];
            res |= (mask & (wrd << shiftLeft));
        }
        return res;
    }

    bool insertValIntoHeap(const std::vector<uint64_t> &dlta,
                           uint64_t startIdx,
                           uint64_t *vecPtr,
                           uint64_t width) {
        uint64_t startBit{0};
        for (uint64_t i = startIdx; i < dlta.size(); i++) {
            uint64_t val = dlta[i];
            if (width == slotWidth and val >= numSamples) {
                std::string msg = "delta index is larger than num_samples. val:"+
                        std::to_string(val)+
                                  " num_samples:" +std::to_string(numSamples);
                throw DeltaManagerException(msg);
            }

            uint64_t &wrd = vecPtr[startBit / 64];
            uint64_t startBitInwrd = startBit % 64;
            wrd |= (val << startBitInwrd);
            uint64_t shiftRight = 64 - startBitInwrd;

            if (shiftRight < width) {
                auto &nextWrd = vecPtr[startBit / 64 + 1];
                nextWrd |= (val >> shiftRight);
            }
            startBit += width;
        }
        return true;
    }

    bool getValFromHeap(std::vector<uint64_t> &dlta,
                        uint64_t cnt,
                        uint64_t *vecPtr,
                        uint64_t width) {
        uint64_t startBit{0};
        uint64_t mask = width < 64? (((uint64_t)1 << width) - 1) : -1;
        for (uint64_t i = 0; i < cnt; i++) {
            uint64_t res{0};
            uint64_t wrd = vecPtr[startBit / 64];
            uint64_t startBitInwrd = startBit % 64;
            res |= (mask & (wrd >> startBitInwrd));
            uint64_t shiftRight = 64 - startBitInwrd;
            if (shiftRight < width) {
                wrd = vecPtr[startBit / 64 + 1];
                res |= (mask & (wrd << shiftRight));
            }
            startBit += width;
            dlta.push_back(res);
        }
        return true;
    }

};

#endif //MANTIS_DELTAMANAGER_H
