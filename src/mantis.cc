/*
 * ============================================================================
 *       Filename:  mantis.cc
 * ============================================================================
 */

#include <iostream>
#include <string>
#include <vector>
#include <cassert>
#include <exception>

#include "MantisFS.h"
#include "ProgOpts.h"
#include "clipp.h"
#include "spdlog/spdlog.h"
#include "mantisconfig.hpp"
#include "mst.h"
#include "mstQuery.h"
#include <omp.h>

#include "tbb/global_control.h"

template <typename T>
void explore_options_verbose(T& res) {
  if(res.any_error()) { std::cerr << "error\n"; }

  //aggregated errors
  if(res.unmapped_args_count()) { std::cerr << "error unmapped args count\n"; /* ... */ }
  if(res.any_bad_repeat()) { std::cerr << "error bad repeat \n"; /* ... */ }
  if(res.any_blocked())    { std::cerr << "error blocked \n"; /* ... */ }
  if(res.any_conflict())   { std::cerr << "error conflict\n"; /* ... */ }

  for(const auto& m : res.missing()) { 
    std::cerr << "missing " << m.param() << " after index " << m.after_index() << '\n';
  }

  //per-argument mapping
  for(const auto& m : res) {
    std::cerr << m.index() << ": " << m.arg() << " -> " << m.param();
    std::cerr << " repeat #" << m.repeat();
    if(m.blocked()) std::cerr << " blocked";
    if(m.conflict()) std::cerr << " conflict";
    std::cerr<< '\n';
  }
}

//int query_main (QueryOpts& opt);
int build_main (BuildOpts& opt);
int build_blockedCQF_main (BuildOpts& opt);
int validate_main (ValidateOpts& opt);
int build_mst_main (QueryOpts& opt);
int mst_query_main(QueryOpts &opt);
int query_main (QueryOpts& opt);
int validate_mst_main(MSTValidateOpts &opt);
int stats_main(StatsOpts& statsOpts);

int merge_main(MergeOpts &opt);
int construct_mantis_by_merge_main(BuildOpts& opt);
int validate_merge_main(ValidateMergeOpts &opt);
int compare_indices_main(CompareIndicesOpt &opt);
int lsmt_initialize_main(LSMT_InitializeOpts &opt);
int lsmt_update_main(LSMT_UpdateOpts &opt);
int lsmt_query_main(LSMT_QueryOpts &opt);



/*
 * ===  FUNCTION  =============================================================
 *         Name:  main
 *  Description:
 * ============================================================================
 */
int main ( int argc, char *argv[] ) {
  using namespace clipp;
  enum class mode {build, build_mst, build_by_merge, validate_mst, query, validate, stats, merge, validate_merge,
                  compare_indices, lsmt_init, lsmt_update, lsmt_query, help};
  mode selected = mode::help;

  auto console = spdlog::stdout_color_mt("mantis_console");

  BuildOpts bopt;
  QueryOpts qopt;
  ValidateOpts vopt;
  MSTValidateOpts mvopt;
  StatsOpts sopt;
  MergeOpts mopt;
  ValidateMergeOpts vmopt;
  CompareIndicesOpt ciopt;
  LSMT_InitializeOpts lsmtiopt;
  LSMT_UpdateOpts lsmtuopt;
  LSMT_QueryOpts lsmtqopt;
  bopt.console = console;
  qopt.console = console;
  vopt.console = console;
  mvopt.console = console;
  sopt.console = console;
  mopt.console = console;
  vmopt.console = console;
  ciopt.console = console;
  lsmtiopt.console = console;
  lsmtuopt.console = console;
  lsmtqopt.console = console;


  auto ensure_file_exists = [](const std::string& s) -> bool {
    bool exists = mantis::fs::FileExists(s.c_str());
    if (!exists) {
      std::string e = "The required input file " + s + " does not seem to exist.";
      throw std::runtime_error{e};
    }
    return true;
  };

  auto ensure_dir_exists = [](const std::string& s) -> bool {
    bool exists = mantis::fs::DirExists(s.c_str());
    if (!exists) {
      std::string e = "The required input directory " + s + " does not seem to exist.";
      throw std::runtime_error{e};
    }
    return true;
  };

  auto build_mode = (
                     command("build").set(selected, mode::build),
                     option("-e", "--eqclass_dist").set(bopt.flush_eqclass_dist) % "write the eqclass abundance distribution",
										 required("-s","--log-slots") & value("log-slots",
																											 bopt.qbits) % "log of number of slots in the output CQF",
                     option("-t", "--threads") & value("num_threads", bopt.numthreads) % "number of threads",
                     required("-i", "--input-list") & value(ensure_file_exists, "input_list", bopt.inlist) % "file containing list of input filters",
                     required("-o", "--output") & value("build_output", bopt.out) % "directory where results should be written"
                     );
  auto build_by_merge_mode = (
          command("build_by_merge").set(selected, mode::build_by_merge),
                  option("-e", "--eqclass_dist").set(bopt.flush_eqclass_dist) % "write the eqclass abundance distribution",
                  required("-s","--log-slots") & value("log-slots",
                                                       bopt.qbits) % "log of number of slots in the output CQF",
                  option("-t", "--threads") & value("num_threads", bopt.numthreads) % "number of threads",
                  option("-p", "--processes") & value("num_processes", bopt.numProcesses) % "number of linux commands to be run at the same time",
                  required("-i", "--input-list") & value(ensure_file_exists, "input_list", bopt.inlist) % "file containing list of input filters",
                  required("-o", "--output") & value("build_output", bopt.out) % "directory where results should be written"
  );

  auto build_mst_mode = (
          command("mst").set(selected, mode::build_mst),
                  required("-p", "--index-prefix") & value(ensure_dir_exists, "index_prefix", qopt.prefix) % "The directory where the index is stored.",
                  option("-t", "--threads") & value("num_threads", qopt.numThreads) % "number of threads",
                  (
                          required("-k", "--keep-RRR").set(qopt.keep_colorclasses) % "Keep the previous color class RRR representation."
                          |
                          required("-d", "--delete-RRR").set(qopt.remove_colorClasses) % "Remove the previous color class RRR representation."
                  )
                  );

  auto validate_mst_mode = (
                  command("validatemst").set(selected, mode::validate_mst),
                  required("-p", "--index-prefix") & value(ensure_dir_exists, "index_prefix", mvopt.prefix) % "The directory where the index is stored.",
                  required("-n", "--num-experiments") & value("num_experiments", mvopt.numSamples) % "Number of experiments."
  );

  auto query_mode = (
                     command("query").set(selected, mode::query),
                     option("-b", "--bulk").set(qopt.process_in_bulk) % "Process the whole input query file as a bulk.",
                     option("-1", "--use-colorclasses").set(qopt.use_colorclasses)
                     % "Use color classes as the color info representation instead of MST",
                     option("-j", "--json").set(qopt.use_json) % "Write the output in JSON format",
                     option("-k", "--kmer") & value("kmer", qopt.k) % "size of k for kmer.",
                     required("-p", "--input-prefix") & value(ensure_dir_exists, "query_prefix", qopt.prefix) % "Prefix of input files.",
                     option("-o", "--output") & value("output_file", qopt.output) % "Where to write query output.",
                     value(ensure_file_exists, "query", qopt.query_file) % "Prefix of input files."
                     );

  auto validate_mode = (
                     command("validate").set(selected, mode::validate),
                     required("-i", "--input-list") & value(ensure_file_exists, "input_list", vopt.inlist) % "file containing list of input filters",
                     required("-p", "--input-prefix") & value(ensure_dir_exists, "dbg_prefix", vopt.prefix) % "Directory containing the mantis dbg.",
                     value(ensure_file_exists, "query", vopt.query_file) % "Query file."
                     );
    auto stats_mode = (
            command("stats").set(selected, mode::stats),
                    required("-p", "--index-prefix") & value(ensure_dir_exists, "index_prefix", sopt.prefix) % "The directory where the index is stored.",
                    required("-n", "--num-samples") & value("number_of_samples", sopt.numSamples) % "Number of experiments.",
                    option("-t", "--type") & value("type", sopt.type) % "what stats? (mono, cc_density, color_dist, jmerkmer), default: mono",
                    option("-j", "--jmer-length") & value("size-of-jmer", sopt.j) % "value of j for constituent jmers of a kmer (default: 23)."
    );

  
  auto merge_mode = (
                    command("merge").set(selected, mode::merge),
										option("-t", "--thread-count") & value("thread-count", mopt.numThreads) % "number of threads to use in intermediate unique color-id filtering phase",
                    required("-i1", "--input-dir-1") & value("input-dir-1", mopt.dir1) % "directory containing the first CdBG",
                    required("-i2", "--input-dir-2") & value("input-dir-2", mopt.dir2) % "directory containing the second CdBG",
                    required("-o", "--output") & value("merge-output", mopt.out) % "directory where the merged CdBG should be written",
                    option("-l", "--time-log").set(mopt.timeLog) % "write the summary time-log for the algorithm",
                    option("-rm", "--remove-indices").set(mopt.removeIndices) % "remove the input mantis indices from disk"
                    );

  auto validate_merge_mode = (
                              command("validate_merge").set(selected, mode::validate_merge),
                              required("-c", "--correct-cdbg") & value("correct-cdbg", vmopt.correctRes) % "directory containing the correct CdBG",
                              required("-m", "--merged-cdbg") & value("merged-cdbg", vmopt.mergeRes) % "directory containing the merged CdBG"
                            );

  auto compare_indices_mode = (
                              command("compare_indices").set(selected, mode::compare_indices),
                              required("-c1", "--cdbg-1") & value("cdbg-1", ciopt.cdbg1) % "directory containing the first CdBG",
                              required("-c2", "--cdbg-2") & value("cdbg-2", ciopt.cdbg2) % "directory containing the second CdBG"
                              );

  auto lsmt_init_mode = (
                        command("lsmt_init").set(selected, mode::lsmt_init),
                        required("-d", "--dir") & value("dir", lsmtiopt.dir)
                        % "directory where the LSM-tree will reside",
                        option("-c", "--scaling-factor") & value("scaling-factor", lsmtiopt.scaling_factor)
                        % "scaling factor for the LSM-tree levels",
                        option("-k", "--kmer-threshold") & value("kmer-threshold", lsmtiopt.kmer_threshold)
                        % "kmer count threshold for the level 0 of the LSM-tree",
                        option("-th", "--cqf-count-threshold") & value("cqf-count-threshold", lsmtiopt.cqf_count_threshold)
                        % "CQF count threshold for the level 0 of the LSM-tree",
                        option("-s", "--sample-threshold") & value("sample-threshold", lsmtiopt.sample_threshold)
                        % "threshold on the count of samples kept pending before insertion into the LSM-tree",
                        option("-q", "--q-bit-init-build") & value("q-bit-init-build", lsmtiopt.qBit_init_build)
                        % "q (quotient)-bit count for the initial mantis build at an update"
                        );

  auto lsmt_update_mode = (
                          command("lsmt_update").set(selected, mode::lsmt_update),
                          required("-d", "--dir") & value("dir", lsmtuopt.dir)
                          % "directory where the LSM-tree resides",
                          required("-i", "--input-list") & value(ensure_file_exists, "input-list", lsmtuopt.input_list)
                          % "file containing list of input sample-filters",
                          option("-t", "--thread-count") & value("thread-count", lsmtuopt.thread_count)
                          % "number of threads to use in intermediate merge operations"
                          );

  auto lsmt_query_mode = (
                          command("lsmt_query").set(selected, mode::lsmt_query),
                          required("-d", "--dir") & value("dir", lsmtqopt.dir)
                          % "directory where the LSM-tree resides",
                          required("-q", "--query-file") & value(ensure_file_exists, "query-file", lsmtqopt.queryFile)
                          % "file containing the query transcripts",
                          required("-o", "--output") & value("query-output", lsmtqopt.output)
                          % "file containing the query results",
                          option("-k", "--kmer-length") & value("kmer-length", lsmtqopt.k) % "length of k-mers"
                        );


  auto cli = (
              (build_mode | build_mst_mode | build_by_merge_mode | validate_mst_mode | query_mode | validate_mode | stats_mode |
              merge_mode | validate_merge_mode | compare_indices_mode |
              lsmt_init_mode | lsmt_update_mode | lsmt_query_mode |
              command("help").set(selected,mode::help) |
               option("-v", "--version").call([]{std::cout << "mantis " << mantis::version << '\n'; std::exit(0);}).doc("show version")
              )
             );

  assert(build_mode.flags_are_prefix_free());
  assert(build_by_merge_mode.flags_are_prefix_free());
  assert(query_mode.flags_are_prefix_free());
  assert(validate_mode.flags_are_prefix_free());
  assert(build_mst_mode.flags_are_prefix_free());
  assert(validate_mst_mode.flags_are_prefix_free());
  assert(stats_mode.flags_are_prefix_free());
  assert(merge_mode.flags_are_prefix_free());
  assert(validate_merge_mode.flags_are_prefix_free());
  assert(compare_indices_mode.flags_are_prefix_free());
  assert(lsmt_init_mode.flags_are_prefix_free());
  assert(lsmt_update_mode.flags_are_prefix_free());
  assert(lsmt_query_mode.flags_are_prefix_free());

  decltype(parse(argc, argv, cli)) res;
  try {
    res = parse(argc, argv, cli);
  } catch (std::exception& e) {
    std::cout << "\n\nParsing command line failed with exception: " << e.what() << "\n";
    std::cout << "\n\n";
    std::cout << make_man_page(cli, "mantis");
    return 1;
  }

  //explore_options_verbose(res);
  int numThreads = 16;
  if (res) {
    switch(selected) {
      case mode::build:numThreads=qopt.numThreads;break;
      case mode::build_by_merge:numThreads=bopt.numthreads;break;
      case mode::merge:numThreads=mopt.numThreads;break;
      case mode::query:numThreads=qopt.numThreads;break;
    }
  }
  tbb::global_control c(tbb::global_control::max_allowed_parallelism, numThreads);
  if(res) {
    switch(selected) {
    case mode::build:
      build_blockedCQF_main(bopt);
      qopt.prefix = bopt.out; qopt.numThreads = bopt.numthreads; qopt.remove_colorClasses = true;
            build_mst_main(qopt);
      break;//build_main(bopt);  break;
    case mode::build_by_merge:
      //        omp_set_dynamic(false);
//        omp_set_num_threads(bopt.numthreads);
        construct_mantis_by_merge_main(bopt);
        break;
    case mode::build_mst: build_mst_main(qopt); break;
    case mode::validate_mst: validate_mst_main(mvopt); break;
    case mode::query: qopt.use_colorclasses? query_main(qopt):mst_query_main(qopt);  break;
    case mode::validate: validate_main(vopt);  break;
    case mode::stats: stats_main(sopt);  break;
    case mode::merge: merge_main(mopt);  break;
    case mode::validate_merge: validate_merge_main(vmopt); break;
    // case mode::compare_indices: compare_indices_main(ciopt); break;
    case mode::lsmt_init: lsmt_initialize_main(lsmtiopt); break;
    case mode::lsmt_update: lsmt_update_main(lsmtuopt); break;
    case mode::lsmt_query: lsmt_query_main(lsmtqopt); break;
    case mode::help: std::cout << make_man_page(cli, "mantis"); break;
    }
  } else {
    auto b = res.begin();
    auto e = res.end();
    if (std::distance(b,e) > 0) {
      if (b->arg() == "build") {
        std::cout << make_man_page(build_mode, "mantis");
      } else if (b->arg() == "mst") {
        std::cout << make_man_page(build_mst_mode, "mantis");
      } else if (b->arg() == "query") {
        std::cout << make_man_page(query_mode, "mantis");
      } else if (b->arg() == "validatemst") {
        std::cout << make_man_page(validate_mst_mode, "mantis");
      } else if (b->arg() == "validate") {
        std::cout << make_man_page(validate_mode, "mantis");
      } else if (b->arg() == "stats") {
        std::cout << make_man_page(stats_mode, "mantis");
      }
        else if (b->arg() == "merge") {
        std::cout << make_man_page(merge_mode, "mantis");
      } else if (b->arg() == "validate_merge") {
        std::cout << make_man_page(validate_merge_mode, "mantis");
      } else if (b->arg() == "compare_indices") {
        std::cout << make_man_page(compare_indices_mode, "mantis");
      } else if (b->arg() == "lsmt_init") {
        std::cout << make_man_page(lsmt_init_mode, "mantis");
      } else if(b->arg() == "lsmt_query") {
        std::cout << make_man_page(lsmt_query_mode, "mantis");
      }
      else {
        std::cout << "There is no command \"" << b->arg() << "\"\n";
        std::cout << usage_lines(cli, "mantis") << '\n';
      }
    } else {
      std::cout << usage_lines(cli, "mantis") << '\n';
    }
  }
  
  return 0;
}