#ifndef __MANTIS_PROG_OPTS__
#define __MANTIS_PROG_OPTS__
#include <memory>
#include "spdlog/spdlog.h"

class BuildOpts {
 public:
	bool flush_eqclass_dist{false};
	int qbits;
  std::string inlist;
  std::string out;
	int numthreads{1};
  std::shared_ptr<spdlog::logger> console{nullptr};
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


#endif //__MANTIS_PROG_OPTS__
