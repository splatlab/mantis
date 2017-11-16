#ifndef __MANTIS_PROG_OPTS__
#define __MANTIS_PROG_OPTS__

class BuildOpts {
 public:
  std::string inlist;
  std::string cutoffs;
  std::string out;
};

class QueryOpts {
 public:
  std::string prefix;
  std::string output{"samples.output"};
  std::string query_file;
  bool use_json{false};
};

class ValidateOpts {
 public:
  std::string inlist;
  std::string cutoffs;
  std::string prefix;
  std::string query_file;
};


#endif //__MANTIS_PROG_OPTS__
