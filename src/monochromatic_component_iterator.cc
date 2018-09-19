//
// Created by Fatemeh Almodaresi on 6/4/18.
//

#include "monochromatic_component_iterator.h"

#define EQS_PER_SLOT 20000000

uint64_t start_time;

namespace dna {

    /////////////// bases /////////////////
    base operator-(base b) {
        return (base) ((~((uint64_t) b)) & 0x3ULL);
    }

    const base bases[4] = {C, A, T, G};
    const std::map<char, base> base_from_char = {{'A', A},
                                                 {'C', C},
                                                 {'G', G},
                                                 {'T', T},
                                                 {'N', A}};
    const std::map<base, char> base_to_char = {{A, 'A'},
                                               {C, 'C'},
                                               {G, 'G'},
                                               {T, 'T'}};

    ///////////// kmers /////////////////////
    kmer::kmer(void) : len(0), val(0) {}

    kmer::kmer(base b) : len(1), val((uint64_t) b) {}

    kmer::kmer(int l, uint64_t v) : len(l), val(v & BITMASK(2 * l)) {
        assert(l <= 32);
    }

    static uint64_t string_to_kmer_val(std::string s) {
        uint64_t val = 0;
        for (auto c : s)
            val = (val << 2) | ((uint64_t) (base_from_char.at(c)));
        return val;
    }

    kmer::kmer(std::string s) : len(s.size()), val(string_to_kmer_val(s)) {
        assert(s.size() <= 32);
    }

    // Convert to string
    kmer::operator std::string() const {
        std::string s;
        for (auto i = 1; i < len + 1; i++)
            s = s + base_to_char.at((base) ((val >> (2 * (len - i))) & BITMASK(2)));
        return s;
    }

    bool operator<(kmer a, kmer b) {
        return a.len != b.len ? a.len < b.len : a.val < b.val;
    }

    bool operator==(kmer a, kmer b) {
        return a.len == b.len && a.val == b.val;
    }

    bool operator!=(kmer a, kmer b) {
        return !operator==(a, b);
    }

    // Return the reverse complement of k
    kmer operator-(kmer k) {
        uint64_t val = k.val;
        val =
                (val >> 32) |
                (val << 32);
        val =
                ((val >> 16) & 0x0000ffff0000ffff) |
                ((val << 16) & 0xffff0000ffff0000);
        val =
                ((val >> 8) & 0x00ff00ff00ff00ff) |
                ((val << 8) & 0xff00ff00ff00ff00);
        val =
                ((val >> 4) & 0x0f0f0f0f0f0f0f0f) |
                ((val << 4) & 0xf0f0f0f0f0f0f0f0);
        val =
                ((val >> 2) & 0x3333333333333333) |
                ((val << 2) & 0xcccccccccccccccc);
        val = ~val;
        val >>= 64 - 2 * k.len;
        return kmer(k.len, val);
    }

    // backwards from standard definition to match kmer.h definition
    kmer canonicalize(kmer k) {
        return -k < k ? k : -k;
    }

    // Return the kmer of length |a| that results from shifting b into a
    // from the right
    kmer operator<<(kmer a, kmer b) {
        uint64_t val = ((a.val << (2 * b.len)) | b.val) & BITMASK(2 * a.len);
        return kmer(a.len, val);
    }

    // Return the kmer of length |b| that results from shifting b into a
    // from the left
    kmer operator>>(kmer a, kmer b) {
        uint64_t val
                = ((b.val >> (2 * a.len)) | (a.val << (2 * (b.len - a.len))))
                  & BITMASK(2 * b.len);
        return kmer(b.len, val);
    }

    // Append two kmers
    kmer operator+(kmer a, kmer b) {
        int len = a.len + b.len;
        assert(len <= 32);
        uint64_t val = (a.val << (2 * b.len)) | b.val;
        return kmer(len, val);
    }

    kmer prefix(kmer k, int len) { return kmer(len, k.val >> (2 * (k.len - len))); }

    kmer suffix(kmer k, int len) { return kmer(len, k.val & BITMASK(2 * len)); }

    bool period_divides(kmer k, uint64_t periodicity) {
        static const uint64_t multipliers[33] =
                {
                        0,
                        0x5555555555555555, // 1
                        0x1111111111111111, // 2
                        0x1041041041041041, // 3
                        0x0101010101010101, // 4
                        0x1004010040100401, // 5
                        0x1001001001001001, // 6
                        0x0100040010004001, // 7
                        0x0001000100010001, // 8
                        0x0040001000040001, // 9
                        0x1000010000100001, // 10
                        0x0000100000400001, // 11
                        0x0001000001000001, // 12
                        0x0010000004000001, // 13
                        0x0100000010000001, // 14
                        0x1000000040000001, // 15
                        0x0000000100000001, // 16
                        0x0000000400000001, // 17
                        0x0000001000000001, // 18
                        0x0000004000000001, // 19
                        0x0000010000000001, // 20
                        0x0000040000000001, // 21
                        0x0000100000000001, // 22
                        0x0000400000000001, // 23
                        0x0001000000000001, // 24
                        0x0004000000000001, // 25
                        0x0010000000000001, // 26
                        0x0040000000000001, // 27
                        0x0100000000000001, // 28
                        0x0400000000000001, // 29
                        0x1000000000000001, // 30
                        0x4000000000000001, // 31
                        0x0000000000000001, // 32
                };
        uint64_t piece = k.val & BITMASK(2 * periodicity);
        piece = piece * multipliers[periodicity];
        piece = piece & BITMASK(2 * k.len);
        return piece == k.val;
    }

    uint64_t period(kmer k) {
        for (int i = 1; i <= k.len; i++) {
            if (period_divides(k, i))
                return i;
        }
        abort();
    }

    canonical_kmer::canonical_kmer(void) : kmer() {}

    canonical_kmer::canonical_kmer(base b) : kmer(canonicalize(kmer(b))) {}

    canonical_kmer::canonical_kmer(int l, uint64_t v)
            : kmer(canonicalize(kmer(l, v))) {}

    canonical_kmer::canonical_kmer(std::string s) : kmer(canonicalize(kmer(s))) {}

    canonical_kmer::canonical_kmer(kmer k) : kmer(canonicalize(k)) {}

}


//////////////////////////////// monochromatic_component_iterator //////////////////////////

monochromatic_component_iterator::monochromatic_component_iterator(const CQF<KeyObject> *g,
                                                                   std::vector<sdsl::rrr_vector < 63>>

&bvin,
uint64_t num_samplesin
)
:

cqf (g), it(g->begin(0)), bvs(bvin), num_samples(num_samplesin) {
    // initialize cqf iterator
    k = cqf->keybits() / 2; // 2-bit encoded
    std::cerr << "k : " << k << "\n";
    sdsl::util::assign(visited, sdsl::bit_vector(cqf->capacity(), 0));
    std::cerr << "kmers: " << cqf->size() << "\n";
    std::cerr << "slots: " << cqf->capacity() << "\n";
    withMax0.resize(9);
    sdsl::bit_vector d(num_samples, 0);
    eqclass_map[d] = 0;
}

monochromatic_component_iterator::work_item
monochromatic_component_iterator::front(std::queue<monochromatic_component_iterator::work_item> &w) {
    return w.front();
}

bool monochromatic_component_iterator::done() { return it.done(); }

void monochromatic_component_iterator::operator++(void) {

    if (it.done()) return; // don't cross the bound (undefined behaviour)
    ++it;
    auto keyFromKmer = [this](KeyObject keyobj) {
        return HashUtil::hash_64i(keyobj.key, BITMASK(this->cqf->keybits()));
    };
    while (!it.done()) {
        if (visitedKeys.find(keyFromKmer(*it)) != visitedKeys.end()) {
            //if ((bool)(visited[it.iter.current])) {
            ++it;
        } else {
            break;
        }
    }
}

Mc_stats monochromatic_component_iterator::operator*(void) {
    if (!work.empty()) {
        std::cerr << "Throw Exception. The work queue should be empty at this point.\n";
        std::exit(1);
    }
    Mc_stats res;
    if (it.done()) return res;

    KeyObject keyobj = *it;
    node root(k, HashUtil::hash_64i(keyobj.key, BITMASK(cqf->keybits())));
    monochromatic_component_iterator::work_item neww = {root, it.iter.current, keyobj.count};
    work.push(neww);

    while (!work.empty()) {
        work_item w = front(work);
        work.pop();
        // pass over those that have been already visited
        if (visitedKeys.find(w.curr.val) != visitedKeys.end()) {
            //if (visited[w.idx]) {
            continue;
        }
        //std::cerr << "for w " << std::string(w.curr) << " : ";
        for (auto &neighbor : neighbors(w)) {
            //std::cerr << "n " << std::string(neighbor.curr) << " ";
            if (neighbor.curr != w.curr) {
                /*if (neighbor.colorid != w.colorid) {
                    res.min_dist = std::min(res.min_dist, manhattanDist(neighbor.colorid, w.colorid));
                }*/
                if (visitedKeys.find(neighbor.curr.val) == visitedKeys.end()) {
                    //if (visited[neighbor.idx] == 0) {
                    if (neighbor.colorid == w.colorid) {
                        work.push(neighbor);
                    }
                }
            }
        }
        //std::cerr << "\n";
        visitedKeys.insert(w.curr.val);
        res.nodeCnt++;
        visited[w.idx] = 1; // set the corresponding bit
        cntr++;
        if (visitedKeys.size() % 1000000 == 0)
            //if (cntr % 1000000 == 0)
            std::cerr << "visited " << cntr << " kmers\n";
        //std::cerr << " idx " << w.idx << " visited " << visited[w.idx] << "\n";

    }
    return res;
}

std::set<monochromatic_component_iterator::work_item>
monochromatic_component_iterator::neighbors(monochromatic_component_iterator::work_item n) {
    std::set<work_item> result;
    for (const auto b : dna::bases) {
        uint64_t eqid, idx;
        if (exists(b >> n.curr/*b + n.curr*/, idx, eqid))
            if (eqid != n.colorid) // ignore the neighbor if it's a self-loop
                result.insert(work_item(b >> n.curr, idx, eqid));
        if (exists(n.curr << b/*n.curr + b*/, idx, eqid))
            if (eqid != n.colorid)
                result.insert(work_item(n.curr << b, idx, eqid));
    }
    return result;
}

bool monochromatic_component_iterator::exists(edge e, uint64_t &idx, uint64_t &eqid) {
    uint64_t tmp = e.val;
    KeyObject key(HashUtil::hash_64(tmp, BITMASK(cqf->keybits())), 0, 0);
    auto eq_idx = cqf->queryValAndIdx(key);
    if (eq_idx.first) {
        eqid = eq_idx.first - 1;
        idx = eq_idx.second;
        return true;
    }
    return false;
}

void monochromatic_component_iterator::buildColor(std::vector<uint64_t> &eq, uint64_t eqid) {
    uint64_t i{0}, bitcnt{0}, wrdcnt{0};
    uint64_t idx = eqid / EQS_PER_SLOT;
    uint64_t offset = eqid % EQS_PER_SLOT;
    //std::cerr << eqid << " " << num_samples << " " << idx << " " << offset << "\n";
    while (i < num_samples) {
        bitcnt = std::min(num_samples - i, (uint64_t) 64);
        uint64_t wrd = (bvs[idx]).get_int(offset * num_samples + i, bitcnt);
        eq[wrdcnt++] = wrd;
        i += bitcnt;
    }
}

uint64_t monochromatic_component_iterator::manhattanDist(uint64_t eqid1, uint64_t eqid2) {
    uint64_t dist{0};
    std::vector<uint64_t> eq1(((num_samples - 1) / 64) + 1), eq2(((num_samples - 1) / 64) + 1);
    buildColor(eq1, eqid1);
    buildColor(eq2, eqid2);

    for (uint64_t i = 0; i < eq1.size(); i++) {
        if (eq1[i] != eq2[i])
            dist += sdsl::bits::cnt(eq1[i] ^ eq2[i]);
    }
    return dist;

}

__uint128_t monochromatic_component_iterator::manhattanDistBvHash(uint64_t eqid1,
                                                                  uint64_t eqid2,
                                                                  uint64_t num_samples = 2586) {
    sdsl::bit_vector dist(num_samples, 0);
    std::vector<uint64_t> eq1(((num_samples - 1) / 64) + 1), eq2(((num_samples - 1) / 64) + 1);
    buildColor(eq1, eqid1);
    buildColor(eq2, eqid2);

    for (uint64_t i = 0; i < eq1.size(); i++) {
        uint64_t bitcnt = std::min(this->num_samples - (i * 64), (uint64_t) 64);
        dist.set_int((i * 64), (eq1[i] ^ eq2[i]), bitcnt);
        //std::cerr << i << " " << eq1[i] << " " << eq2[i] << " " << (eq1[i] ^ eq2[i]) << "\n";
    }
    __uint128_t dist_hash = HashUtil::MurmurHash128A((void *) dist.data(),
                                                     dist.capacity() / 8, 2038074743,
                                                     2038074751);
    return dist_hash;
}

void monochromatic_component_iterator::manhattanDistBvHash(uint64_t eqid1,
                                                           uint64_t eqid2,
                                                           sdsl::bit_vector &dist,
                                                           uint64_t num_samples = 2586) {
    std::vector<uint64_t> eq1(((num_samples - 1) / 64) + 1), eq2(((num_samples - 1) / 64) + 1);
    buildColor(eq1, eqid1);
    buildColor(eq2, eqid2);

    for (uint64_t i = 0; i < eq1.size(); i++) {
        uint64_t bitcnt = std::min(this->num_samples - (i * 64), (uint64_t) 64);
        dist.set_int((i * 64), (eq1[i] ^ eq2[i]), bitcnt);
        //std::cerr << i << " " << eq1[i] << " " << eq2[i] << " " << (eq1[i] ^ eq2[i]) << "\n";
    }
}


void monochromatic_component_iterator::neighborDist(uint64_t cntrr) {
    KeyObject keyobj = *it;
/*
    if (keyobj.count == 0)
        std::cerr << "we have 0\n";
    if (keyobj.count >= 19216547)
        std::cerr << cntrr << ", keyobj cnt: " << keyobj.count << "\n";
*/
    node curn(k, HashUtil::hash_64i(keyobj.key, BITMASK(cqf->keybits())));
    work_item cur = {curn, it.iter.current, keyobj.count - 1};
    uint64_t mind{UINTMAX_MAX}, meand{0}, maxd{0}, neighborCnt{0};
    for (auto &nei : neighbors(cur)) {
        neighborCnt++;
        if (nei.colorid != cur.colorid) {
            auto d = manhattanDist(nei.colorid, cur.colorid);
            mind = std::min(mind, d);
            maxd = std::max(maxd, d);
            meand += d;
        } else {
            mind = 0;
        }
    }
    // if the node is isolated (or only has a self-loop) it should have maximum possible distance
    if (neighborCnt == 0) {
        isolatedCnt++;
        return;
    }

    if (!maxd) {
        withMax0[neighborCnt]++;
        return;
    }
    // when we get here, neighborCnt is > 0, and mind is not UINT_MAX
    std::cout << neighborCnt << "\t" << mind << "\t"
              << (meand / neighborCnt) << "\t" << maxd << "\n";
}

void monochromatic_component_iterator::buildEqGraph(uint64_t cntrr) {
    KeyObject keyobj = *it;
    node curn(k, HashUtil::hash_64i(keyobj.key, BITMASK(cqf->keybits())));
    work_item cur = {curn, it.iter.current, keyobj.count - 1};
    uint64_t neighborCnt{0};
    for (auto &nei : neighbors(cur)) {
        neighborCnt++;
        if (nei.colorid < cur.colorid) {
            uint16_t d = (uint16_t) manhattanDist(nei.colorid, cur.colorid);
            //if (d <= distThreshold) {
            Edge e(static_cast<uint32_t>(nei.colorid), static_cast<uint32_t>(cur.colorid));
            if (edges.find(e) == edges.end()) {
                edges[e] = d;
            } else if (edges[e] > d) {
                edges[e] = d;
            }
            //}
        }
    }

}

void monochromatic_component_iterator::uniqNeighborDist(uint64_t num_samples) {
    KeyObject keyobj = *it;
    node curn(k, HashUtil::hash_64i(keyobj.key, BITMASK(cqf->keybits())));
    work_item cur = {curn, it.iter.current, keyobj.count - 1};
    uint64_t mind{UINTMAX_MAX}, meand{0}, maxd{0}, neighborCnt{0};
    //std::cerr << " current : " << cur.colorid << "\n";
    for (auto &nei : neighbors(cur)) {
        neighborCnt++;
        sdsl::bit_vector d(num_samples, 0);
        if (nei.colorid != cur.colorid) {
            //xstd::cerr << nei.colorid << "\n";
            //auto d = manhattanDistBvHash(nei.colorid, cur.colorid, num_samples);
            manhattanDistBvHash(nei.colorid, cur.colorid, d, num_samples);
            if (eqclass_map.find(d) == eqclass_map.end()) {
                eqclass_map[d] = 1;
            } else {
                eqclass_map[d]++;
            }
        } /*else {
            eqclass_map[d]++;
        }*/
    }
}

/*
 * ===  FUNCTION  ============================================================
 *         Name:  main
 *  Description:
 * ===========================================================================
 */

int main(int argc, char *argv[]) {

    std::string command = argv[1];
    std::string cqf_file = argv[2];
    std::string eqlistfile = argv[3];
    uint64_t num_samples = std::stoull(argv[4]);
    //uint64_t num_samples = 2586;
    if (argc > 4)
        num_samples = std::stoull(argv[4]);
    std::cerr << "num samples: " << num_samples << "\n";
    CQF<KeyObject> cqf(cqf_file, true);
    std::cerr << "cqf loaded: " << cqf.size() << "\n";
    std::string eqfile;
    std::ifstream eqlist(eqlistfile);
    std::vector<sdsl::rrr_vector < 63>>
    bvs;
    bvs.reserve(20);
    if (eqlist.is_open()) {
        uint64_t accumTotalEqCls = 0;
        while (getline(eqlist, eqfile)) {
            sdsl::rrr_vector<63> bv;
            bvs.push_back(bv);
            sdsl::load_from_file(bvs.back(), eqfile);
        }
    }
    //BitVectorRRR bv(eqfile);
    std::cerr << "num eq clss: " << ((bvs.size() - 1) * EQS_PER_SLOT * num_samples + bvs.back().size()) / num_samples
              << "\n";
    monochromatic_component_iterator mci(&cqf, bvs, num_samples);
    if (command == "monocomp") {
        while (!mci.done()) {
            std::cout << (*mci).nodeCnt << "\n";
            ++mci;
        }
    } else if (command == "neighborDist") {
        size_t cntrr = 0;
        while (!mci.done()) {
            mci.neighborDist(cntrr++);
            ++mci;
            if (cntrr % 10000000 == 0) {
                std::cerr << cntrr << " done\n";
            }
        }
        std::cerr << "\n\nIn order of 0 to 8:\n";
        uint64_t total0s{0};
        for (auto val : mci.withMax0) {
            std::cerr << val << ",";
            total0s += val;
        }
        std::cerr << "\ntotal 0s: " << total0s - mci.withMax0[0]
                  << "\ntota isolated kmers: " << mci.isolatedCnt
                  << "\n";
    } else if (command == "buildEqGraph") {
        size_t cntrr = 0;
        while (!mci.done()) {
            mci.buildEqGraph(cntrr++);
            ++mci;
            if (cntrr % 10000000 == 0) {
                std::cerr << cntrr << " kmers & " << mci.edges.size() << " edges\n";
            }
        }

        for (auto &kv : mci.edges) {
            std::cout << kv.first.n1 << "\t" << kv.first.n2 << "\t" << kv.second << "\n";
        }
    } else if (command == "uniquDistanceDistribution") {
        uint64_t cntr{0};
        while (!mci.done()) {
            mci.uniqNeighborDist(num_samples);
            if (++cntr % 1000000 == 0) {
                std::cerr << cntr << " kmers & " << mci.eqclass_map.size() << " unique distances.\n";
            }
            ++mci;
        }
        std::cerr << "writing the list of distinct distances:\n";
        sdsl::bit_vector outbv(mci.eqclass_map.size() * num_samples, 0);
        std::cerr << "total cnt: " << mci.eqclass_map.size() << " "
                  << (mci.eqclass_map.size() * num_samples) << "\n";
        cntr = 0;
        for (auto &eq_keyval : mci.eqclass_map) {
            std::cout << eq_keyval.second << "\n";
            auto &bv = eq_keyval.first;
            uint64_t b = 0;
            while (b < num_samples) {
                uint64_t bitcnt = std::min(num_samples - b, (uint64_t) 64);
                uint64_t wrd = bv.get_int(b, bitcnt);
                outbv.set_int(num_samples * cntr + b, wrd, bitcnt);
                b += bitcnt;
            }
            cntr++;
            //std::cerr << cntr << "\n";
        }
        uint64_t till = cqf_file.find_last_of('/');
        sdsl::rrr_vector<> cbv(outbv);
        sdsl::store_to_file(cbv, cqf_file.substr(0, till + 1) + "dist_bv.rrr");
    }
}
