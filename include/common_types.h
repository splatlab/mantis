#ifndef __MANTIS_COMMON_TYPES__
#define __MANTIS_COMMON_TYPES__

#include <unordered_set>
#include <set>
#include <unordered_map>
#include <vector>

namespace mantis {
  using KmerHash = uint64_t;
  using ExperimentID = uint64_t;
  using QuerySet = std::unordered_set<KmerHash>;
//  using QuerySet = std::set<KmerHash>;
  using QuerySets = std::vector<QuerySet>;
  using QueryMap = std::unordered_map<KmerHash, std::pair<uint64_t, std::vector<uint32_t>>>;
  using EqMap = std::unordered_map<KmerHash, std::vector<ExperimentID>>;
  struct BulkQuery {
    QuerySet qs;
     EqMap qmap;
  };


  using QueryResult = std::vector<uint64_t>;//std::unordered_map<uint64_t, uint64_t>;
  using QueryResults = std::vector<QueryResult>;
}

#endif //__MANTIS_COMMON_TYPES__
