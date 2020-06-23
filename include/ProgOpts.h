#ifndef __MANTIS_PROG_OPTS__
#define __MANTIS_PROG_OPTS__
#include <memory>
#include "spdlog/spdlog.h"
#include "json.hpp"
#include "mantisconfig.hpp"


class BuildOpts {
 public:
	bool flush_eqclass_dist{false};
	int qbits;
  std::string inlist;
  std::string out;
  uint32_t numthreads{1};
  uint32_t numProcesses{1};
  std::shared_ptr<spdlog::logger> console{nullptr};

  nlohmann::json to_json() {
    nlohmann::json j;
    j["dump_eqclass_dist"] = flush_eqclass_dist;
    j["quotient_bits"] = qbits;
    j["input_list"] = inlist;
    j["output_dir"] = out;
    j["num_threads"] = numthreads;
    return j;
  }
};

class QueryOpts {
 public:
  std::string prefix;
  std::string output{"samples.output"};
  std::string query_file;
  uint64_t k = 0;
  uint32_t numThreads = 1;
  bool use_json{false};
  std::shared_ptr<spdlog::logger> console{nullptr};
  bool process_in_bulk{false};
  bool use_colorclasses{false};
  bool keep_colorclasses{false};
  bool remove_colorClasses{false};
};

class ValidateOpts {
 public:
  std::string inlist;
  std::string prefix;
  std::string query_file;
  std::shared_ptr<spdlog::logger> console{nullptr};
};

class MSTValidateOpts {
public:
    std::string prefix;
    std::uint64_t numSamples;
    std::uint16_t k;
    std::shared_ptr<spdlog::logger> console{nullptr};
};

class StatsOpts {
public:
    std::string prefix;
    std::string type = "mono";
    std::uint64_t numSamples;
    std::uint64_t j = 23;
    std::shared_ptr<spdlog::logger> console{nullptr};
};



class MergeOpts
{
	public:
		uint numThreads{1};
		bool timeLog{false};
		bool removeIndices{false};
		std::string dir1;
		std::string dir2;
		std::string out;
		std::shared_ptr<spdlog::logger> console{nullptr};

		nlohmann::json to_json() {
			nlohmann::json j;
			j["output_dir"] = out;
			j["input_mantis1"] = dir1;
			j["input_mantis2"] = dir2;
			j["num_threads"] = numThreads;
			j["remove_input_indices"] = removeIndices;
			return j;
		}
};



class ValidateMergeOpts
{
	public:
    std::string correctRes;
		std::string mergeRes;
		std::shared_ptr<spdlog::logger> console{nullptr};
};



class CompareIndicesOpt
{
	public:
    std::string cdbg1;
		std::string cdbg2;
		std::shared_ptr<spdlog::logger> console{nullptr};
};



class LSMT_InitializeOpts
{
  public:
    std::string dir;
    uint scaling_factor{mantis::SCALING_FACTOR};
    uint64_t kmer_threshold{mantis::KMER_THRESHOLD};
    uint32_t cqf_count_threshold{mantis::CQF_COUNT_THRESHOLD};
    uint64_t sample_threshold{mantis::SAMPLE_THRESHOLD};
    uint qBit_init_build{mantis::INIT_BUILD_Q_BIT};
    std::shared_ptr<spdlog::logger> console{nullptr};
};



class LSMT_UpdateOpts
{
  public:
    std::string dir;
    std::string input_list;
    uint thread_count{1};
    std::shared_ptr<spdlog::logger> console{nullptr};
};



class LSMT_QueryOpts
{
  public:
    std::string dir;
    std::string queryFile;
    std::string output;
    uint64_t k = 0;
    // uint32_t numThreads = 1;
    // bool use_json{false};
    bool process_in_bulk{false};
    // bool use_colorclasses{false};
    // bool keep_colorclasses{false};
    // bool remove_colorClasses{false};
    std::shared_ptr<spdlog::logger> console{nullptr};
};

#endif //__MANTIS_PROG_OPTS__