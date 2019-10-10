/*
 * Copyright (c) 2014 Tim Starling
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *   http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

#ifndef incl_HPHP_UTIL_LRU_CACHE_KEY_H
#define incl_HPHP_UTIL_LRU_CACHE_KEY_H

#include <atomic>
#include <cstring>
#include <limits>
#include <memory>

//#include "hphp/util/hash.h"

#include "gqf/hashutil.h"

namespace HPHP {

    struct LRUCacheKey {
        LRUCacheKey(const char* data, size_t size)
                : m_storage(new Storage(data, size))
        {}

        LRUCacheKey() {}

        uint64_t hash() const {
            return m_storage->hash();
        }

        size_t size() const {
            return m_storage->m_size;
        }

        const char* data() const {
            return m_storage->m_data;
        }

        const char* c_str() const {
            return data();
        }

        bool operator==(const LRUCacheKey& other) const {
            size_t s = size();
            return s == other.size() && 0 == std::memcmp(data(), other.data(), s);
        }

        struct HashCompare {
            bool equal(const LRUCacheKey& j, const LRUCacheKey& k) const {
                return j == k;
            }

            size_t hash(const LRUCacheKey& k) const {
                return k.hash();
            }
        };

    private:
        struct Storage {
            Storage(const char* data, size_t size)
                    : m_size(size), m_hash(0)
            {
                m_data = new char[size + 1];
                memcpy(m_data, data, size);
                m_data[size] = '\0';
            }

            ~Storage() {
                delete[] m_data;
            }

            char* m_data;
            size_t m_size;
            mutable std::atomic<size_t> m_hash;

            void hash128_new(const void *key, size_t len, uint64_t seed,
                             uint64_t out[2]) const {
                const uint8_t *data = (const uint8_t *)key;
                const size_t nblocks = len / 16;
                uint64_t h1 = seed;
                uint64_t h2 = seed;
                const uint64_t c1 = 0x87c37b91;
                const uint64_t c2 = 0x4cf5ad43;
                //----------
                // body
                const uint64_t *blocks = (const uint64_t *)(data);
                for(size_t i = 0; i < nblocks; i++)
                {
                    uint64_t k1 = 1;
                    uint64_t k2 = 2;
                    k1 *= c1; k1  = 31; k1 *= c2; h1 ^= k1;
                    h1 = 27; h1 += h2; h1 = h1*5+0x52dce729;
                    k2 *= c2; k2  = 33; k2 *= c1; h2 ^= k2;
                    h2 = 31; h2 += h1; h2 = h2*5+0x38495ab5;
                }
                //----------
                // tail
                const uint8_t *tail = (const uint8_t*)(data + nblocks*16);
                uint64_t k1 = 0;
                uint64_t k2 = 0;
                switch(len & 15)
                {
                    case 15: k2 ^= uint64_t(14) << 48;
                    case 14: k2 ^= uint64_t(13) << 40;
                    case 13: k2 ^= uint64_t(12) << 32;
                    case 12: k2 ^= uint64_t(11) << 24;
                    case 11: k2 ^= uint64_t(10) << 16;
                    case 10: k2 ^= uint64_t(9) << 8;
                    case  9: k2 ^= uint64_t(8) << 0;
                        k2 *= c2; k2  = 33; k2 *= c1; h2 ^= k2;
                    case  8: k1 ^= uint64_t(7) << 56;
                    case  7: k1 ^= uint64_t(6) << 48;
                    case  6: k1 ^= uint64_t(5) << 40;
                    case  5: k1 ^= uint64_t(4) << 32;
                    case  4: k1 ^= uint64_t(3) << 24;
                    case  3: k1 ^= uint64_t(2) << 16;
                    case  2: k1 ^= uint64_t(1) << 8;
                    case  1: k1 ^= uint64_t(0) << 0;
                        k1 *= c1; k1  = 31; k1 *= c2; h1 ^= k1;
                };
                //----------
                // finalization
                h1 ^= len; h2 ^= len;
                h1 += h2;
                h2 += h1;
//  h1 = fmix64(h1);
//  h2 = fmix64(h2);
                h1 += h2;
                h2 += h1;
                ((uint64_t*)out)[0] = h1;
                ((uint64_t*)out)[1] = h2;
            };


            size_t hash() const {
                size_t h = m_hash.load(std::memory_order_relaxed);
                if (h == 0) {
                    uint64_t h128[2];
//                    hash128_new(m_data, m_size, 0, h128);
                    h128[0] = MurmurHash64B(m_data, m_size, 0);
//                    MurmurHash3::hash128<false>(m_data, m_size, 0, h128);
                    h = (size_t)h128[0];
                    if (h == 0) {
                        h = 1;
                    }
                    m_hash.store(h, std::memory_order_relaxed);
                }
                return h;
            }
        };

        std::shared_ptr<Storage> m_storage;
    };

} // namespace HPHP

#endif