#ifndef __MANTIS_COMMON_TYPES__
#define __MANTIS_COMMON_TYPES__

#include <unordered_set>
#include <unordered_map>
#include <vector>

namespace mantis {
  using KmerHash = uint64_t;
  using ExperimentID = uint64_t;
  using QuerySet = std::unordered_set<KmerHash>;
  using QuerySets = std::vector<QuerySet>;
  struct BulkQuery {
    QuerySet qs;
    std::unordered_map<KmerHash, std::vector<ExperimentID>> qmap;
  };


  using QueryResult = std::unordered_map<uint64_t, uint64_t>;
  using QueryResults = std::vector<QueryResult>;
}

#endif //__MANTIS_COMMON_TYPES__
