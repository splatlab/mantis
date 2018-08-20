//
// Created by Fatemeh Almodaresi on 8/17/18.
//

#ifndef MANTIS_DELTAMANAGER_H
#define MANTIS_DELTAMANAGER_H

#include <vector>
#include <math.h>
#include <assert.h>
#include <iostream>
//#include <bitset>
class DeltaManager {
public:

    DeltaManager(uint64_t numSamples,
                 uint64_t approximateColorClsCnt,
                 uint64_t approximateAvgDeltaPerColorCls) {
        slotWidth = ceil(log2(numSamples));
        if (slotWidth * numSamples < 64) {
            slotsPerColorCls = numSamples + 1;
        } else if (slotWidth * approximateAvgDeltaPerColorCls < 64) {
            slotsPerColorCls = (64 / slotWidth + 1) + 1;
        } else {
            slotsPerColorCls = approximateAvgDeltaPerColorCls + 1; // 1 for storing count of deltas per index
        }
        deltas.reserve(approximateColorClsCnt * slotsPerColorCls);
        // assumption: count of slots * their width is greater than 64 bits
        slotsPerColorClsWithPtrs = ((slotsPerColorCls - 1) * slotWidth - 64) / slotWidth + 1;
        colorCnt = 0;
        mask = (((uint64_t) 1) << slotWidth) - 1;
    }

    void insertDeltas(uint64_t colorId, const std::vector<uint64_t> &dlta) {
        // see if we need to split deltas between main DS and heap
        assert(colorId <= colorCnt);
        auto startBit = colorId*slotsPerColorCls*slotWidth;
        auto nextStartBit = (colorId+1)*slotsPerColorCls*slotWidth;
        // We always want to assign slotsPerColorCls slots to each index
        // even if deltas in that index are fewer than the avg num of deltas
        if (colorId == colorCnt) {
            while (deltas.size() < nextStartBit/64+1) {
                deltas.push_back(0);
            }
        }
        // now assuming we've already reserved the space for the new colorId, insert deltas
        uint64_t mainDSDeltaCnt = dlta.size() < slotsPerColorCls ? dlta.size() : slotsPerColorClsWithPtrs - 1;
        insertValIntoDeltaV(startBit, dlta.size(), slotWidth); // insert count of deltas
        startBit += slotWidth;
        for (uint64_t i = 0; i < mainDSDeltaCnt; i++) { // insert values into main DS
            insertValIntoDeltaV(startBit, dlta[i], slotWidth);
            startBit += slotWidth;
        }
        if (mainDSDeltaCnt < dlta.size()) { // in case count of deltas exceeds the reserved count
            // store the rest in an array in heap
            uint64_t *theRest = new uint64_t[(dlta.size() - mainDSDeltaCnt) * slotWidth / 64 + 1];
            insertValIntoDeltaV(startBit, reinterpret_cast<uint64_t>(theRest),
                                64); // store the pointer to the heap in the main DS
            insertValIntoHeap(dlta, mainDSDeltaCnt, theRest, slotWidth);
        }
        if (colorId == colorCnt)
            colorCnt++;
    }

    std::vector<uint64_t> getDeltas(uint64_t colorId) {
        assert(colorId < colorCnt);
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

private:
    uint64_t slotWidth;
    uint64_t slotsPerColorCls;
    uint64_t slotsPerColorClsWithPtrs;
    std::vector<uint64_t> deltas;
    uint64_t colorCnt;
    uint64_t mask;

    // width is limited to 64 (word size)
    bool insertValIntoDeltaV(uint64_t startBit, uint64_t val, uint64_t width) {
        uint64_t localMask = (((uint64_t)1 << width) - 1);
        assert(startBit < deltas.size() * 64);
        uint64_t &wrd = deltas[startBit / 64];
        uint64_t startBitInwrd = startBit % 64;
        wrd &= ~(localMask << startBitInwrd);
        wrd |= (val << startBitInwrd);
        uint64_t shiftRight = 64 - startBitInwrd;
        if (shiftRight < width) {
            assert(startBit < (deltas.size()-1)*64);
            auto &nextWrd = deltas[startBit/64+1];
            nextWrd &= ~(localMask >> shiftRight);
            nextWrd |= (val >> shiftRight);
        }

        return true;
    }

    // width is limited to 64 (word size)
    uint64_t getValFromMDeltaV(uint64_t startBit, uint64_t width) {
        assert(startBit < deltas.size() * 64);
        uint64_t res{0};
        uint64_t wrd = deltas[startBit / 64];
        uint64_t startBitInwrd = startBit % 64;
        res |= (mask & (wrd >> startBitInwrd));
        uint64_t shiftLeft = 64 - startBitInwrd;
        if (shiftLeft < width) {
            assert(startBit < (deltas.size() - 1) * 64);
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
