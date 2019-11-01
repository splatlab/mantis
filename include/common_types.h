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
  using QueryMap = std::unordered_map<KmerHash, uint64_t>;
  using EqMap = std::unordered_map<KmerHash, std::vector<ExperimentID>>;
  struct BulkQuery {
    QuerySet qs;
     EqMap qmap;
  };


  using QueryResult = std::vector<uint64_t>;//std::unordered_map<uint64_t, uint64_t>;
  using QueryResults = std::vector<QueryResult>;


  using Kmer_t = uint64_t;
  using ColorId_t = uint64_t;
  using TranscriptId_t = uint32_t;
  using Abundance_t = uint64_t;
  using SampleId_t = uint64_t;
  // using TranscriptList_t = std::vector<TranscriptId_t>;
  // using KmerTranscriptList_t = std::pair<Kmer_t, TranscriptList_t>;
  // using ColorTranscriptPair_t = std::pair<ColorId_t, TranscriptId_t>;
  // using TranscriptSampleDist_t = std::vector<std::vector<Abundance_t>>;
}

#endif //__MANTIS_COMMON_TYPES__
