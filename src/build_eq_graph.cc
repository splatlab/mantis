#include <iostream>
#include <string>
#include <cmath>
#include <unordered_map>
#include <unordered_set>

#include "bitvector.h"
#include "sdsl/rrr_vector.hpp"
#include "hashutil.h"
#include "clipp.h"

uint32_t hamming_dist(BitVectorRRR& v, uint64_t vec_length, std::vector<uint32_t>& blengths,
                      std::vector<uint32_t>& boffs, uint32_t i, uint32_t j, uint32_t thresh) {
  auto d = [](uint64_t x, uint64_t y) -> int {
      uint64_t res = x ^ y;
      return __builtin_popcountll (res);
  };

  auto idx_i = i * vec_length;
  auto idx_j = j * vec_length;

  uint64_t dist{0};
  uint64_t eqwords{0};
  for (size_t idx = 0; idx < blengths.size(); ++idx) {
    auto wx = v.get_int(idx_i + boffs[idx], blengths[idx]);
    auto wy = v.get_int(idx_j + boffs[idx], blengths[idx]);
    dist += d(wx, wy);
    eqwords += (wx == wy) ? 1 : 0;
    if (dist > thresh) { return dist; }
  }
  if (dist == 0 ) {
      if (eqwords < blengths.size()) {
        std::cerr << "dist = " << dist << ", but eqwords = " << eqwords << " / " << blengths.size() << "\n";
      }
      std::cerr << "for " << i << ", " << j << " eqwords = " << eqwords << "\n";
  }
  return dist;
}

int main(int argc, char* argv[]){
  using namespace clipp; using std::cout; using std::string;
  uint32_t thresh{5};
  std::string fname;//{argv[1]};
  uint64_t vec_length{0};//= std::stoul(argv[2]);
  uint64_t num_to_process{std::numeric_limits<uint64_t>::max()};// = std::stoul(argv[3]);
  bool help{false};
  auto cli = (
              value("input file", fname),
              required("-l", "--length") & value("vec_length", vec_length) % "vector length",
              option("-n") & value("nproc", num_to_process) % "maximum vectors to process (default = all)",
              option("-t") & value("thresh", thresh) % "threshold of distance to report edges (default = 5)",
              option("-h", "--help").set(help)
              );

  if (!parse(argc, argv, cli)) {
    std::cerr << make_man_page(cli, argv[0]);
    std::exit(1);
  } else if (help) {
    std::cerr << make_man_page(cli, argv[0]);
    std::exit(1);
  }

  BitVectorRRR v(fname);

  std::vector<std::unordered_map<uint64_t, std::vector<uint32_t>>> pattern_map;
  auto num_buckets = std::ceil(vec_length/64.0);
  pattern_map.resize(num_buckets);

  std::vector<uint32_t> blengths;
  std::vector<uint32_t> boffs;
  size_t offset{0};
  std::cerr << "[";
  for (size_t i = 0; i < num_buckets; ++i) {
    auto s = std::min(vec_length - offset, static_cast<uint64_t>(64u));
    boffs.push_back(offset);
    blengths.push_back(s);
    std::cerr << " (" << offset << ", " << s << ")";
    offset += s;
  }
  std::cerr << " ]\n";
  std::cerr << "pattern map has " << num_buckets << " buckets\n";

  auto tot_bits = v.bit_size();
  std::cerr << "length = " << tot_bits << "\n";

  offset = 0;
  uint64_t vec_idx{0};
  size_t max_num_vec = num_to_process;
  while (offset < tot_bits and vec_idx < max_num_vec) {
    for (size_t idx = 0; idx < blengths.size(); ++idx) {
      auto w = v.get_int(offset + boffs[idx], blengths[idx]);
      pattern_map[idx][w].push_back(vec_idx);
    }
    ++vec_idx;
    offset = vec_idx * vec_length;
    if ((vec_idx > 0) and (vec_idx % 500000 == 0)) {
      std::cerr << "processed " << vec_idx << " vectors\n";
    }
  }

  uint64_t s{0};
  std::vector<double> means;
  for (size_t b = 0; b < num_buckets; ++b) {
    std::cerr << "bucket " << b << " has " << pattern_map[b].size() << " patterns (" << (static_cast<float>(vec_idx) / pattern_map[b].size()) << ")\n";
    s += std::ceil((static_cast<float>(vec_idx) / pattern_map[b].size()));
    means.push_back((static_cast<float>(vec_idx) / pattern_map[b].size()));
  }

  std::unordered_map<uint64_t, uint32_t> candidates;
  std::unordered_set<uint32_t> skip;
  std::vector<uint64_t> cd;
  cd.reserve(s);
  size_t nn{0};
  size_t nc{0};
  size_t m{std::numeric_limits<size_t>::max()};
  std::cout << max_num_vec << '\t' << vec_length << '\n';
  for (size_t x = 0; x < max_num_vec; ++x) {
    offset = vec_length * x;
    for (size_t idx = 0; idx < blengths.size(); ++idx) {
      auto w = v.get_int(offset + boffs[idx], blengths[idx]);
      auto& c1 = pattern_map[idx][w];
      //cd.insert(cd.end(), c1.begin(), c1.end());
      if (c1.size() < 2*means[idx]) {
        for(auto c : c1) {
          if (c > x) { candidates[c] += 1; }
        }
      }
    }
    //std::sort(cd.begin(), cd.end());
    //auto last = std::unique(cd.begin(), cd.end());
    //nc += std::distance(cd.begin(), last);

    //std::cerr << "x = " << x << ", num candidates = " << candidates.size() << '\n';

    for (auto& kv : candidates) {
      //if (kv.second >= 35) {
        nn++;
        auto d = hamming_dist(v, vec_length, blengths, boffs, x, kv.first, thresh);
        if (d==0) {
          std::cerr << "WHAT?! " << x << " = " << kv.first << ", num matching blocks = " << kv.second << "\n"; 
        for(size_t j = 0; j < vec_length; ++j) {
          std::cerr << v[x*vec_length + j];
        }
        std::cerr << "\n\n";
        for(size_t j = 0; j < vec_length; ++j) {
          std::cerr << v[kv.first*vec_length + j];
        }
        std::cerr << "\n\n";
        }
        bool within_bound = d <= thresh;
        nc += (within_bound);
        m = (d < m) ? d : m;
        if (within_bound) { std::cout << x << '\t' << kv.first << '\t' << d << '\n'; }
        //}
    }
    if (x > 0 and x % 10000 == 0) {
      std::cerr << "x = " << x << " ( " << (nc / static_cast<float>(x)) << ", " << (nn / static_cast<float>(x)) << " : " << m << ")\n";
    }
    candidates.clear();
    //cd.clear();
  }



  return 0;
}
