#ifndef __MANTIS_PROG_OPTS__
#define __MANTIS_PROG_OPTS__
#include <memory>
#include "spdlog/spdlog.h"
#include "json.hpp"


class BuildOpts {
 public:
	bool flush_eqclass_dist{false};
	int qbits;
  std::string inlist;
  std::string out;
	int numthreads{1};
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
  bool use_json{false};
  std::shared_ptr<spdlog::logger> console{nullptr};
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

#endif //__MANTIS_PROG_OPTS__
