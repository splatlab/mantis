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

int query_main(QueryOpts& opt);
int build_dbgccmst(DBGCCMSTOpts& opt);
int build_main(BuildOpts& opt);
int validate_main(ValidateOpts& opt);

int main ( int argc, char *argv[] ) {
  using namespace clipp;
  enum class mode {build, query, dbgccmst, validate, help};
  mode selected = mode::help;

  auto console = spdlog::stdout_color_mt("mantis_console");

  BuildOpts bopt;
  QueryOpts qopt;
  DBGCCMSTOpts mopt;
  ValidateOpts vopt;
  bopt.console = console;
  qopt.console = console;
  mopt.console = console;
  vopt.console = console;

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
                     required("-i", "--input-list") & value(ensure_file_exists, "input_list", bopt.inlist) % "file containing list of input filters",
                     required("-c", "--cutoff-list") & value(ensure_file_exists, "cutoff_list", bopt.cutoffs) % "file containing list of experiment-specific cutoffs",
                     required("-o", "--output") & value("build_output", bopt.out) % "directory where results should be written"
                     );
  auto query_mode = (
                     command("query").set(selected, mode::query),
                     option("-j", "--json").set(qopt.use_json) % "Write the output in JSON format",
                     required("-p", "--input-prefix") & value(ensure_dir_exists, "query_prefix", qopt.prefix) % "Prefix of input files.",
                     option("-o", "--output") & value("output_file", qopt.output) % "Where to write query output.",
                     value(ensure_file_exists, "query", qopt.query_file) % "Prefix of input files."
                     );

  auto dbgccmst_mode = (
                     command("dbgccmst").set(selected, mode::dbgccmst),
                     required("-p", "--input-dir") & value(ensure_dir_exists, "input_dir", mopt.prefix) % "Directory containing input files."
                     );

  auto validate_mode = (
                     command("validate").set(selected, mode::validate),
                     required("-i", "--input-list") & value(ensure_file_exists, "input_list", vopt.inlist) % "file containing list of input filters",
                     required("-c", "--cutoff-list") & value(ensure_file_exists, "cutoff_list", vopt.cutoffs) % "file containing list of experiment-specific cutoffs",
                     required("-p", "--input-prefix") & value(ensure_dir_exists, "dbg_prefix", vopt.prefix) % "Directory containing the mantis dbg.",
                     value(ensure_file_exists, "query", vopt.query_file) % "Query file."
                     );

  auto cli = (
              (build_mode | query_mode | dbgccmst_mode | validate_mode | command("help").set(selected,mode::help) ),
              option("-v", "--version").call([]{std::cout << "version 1.0\n\n";}).doc("show version")  );

  assert(build_mode.flags_are_prefix_free());
  assert(query_mode.flags_are_prefix_free());
  assert(dbgccmst_mode.flags_are_prefix_free());
  assert(validate_mode.flags_are_prefix_free());

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

  if(res) {
    switch(selected) {
    case mode::build: build_main(bopt);  break;
    case mode::query: query_main(qopt);  break;
    case mode::dbgccmst: build_dbgccmst(mopt);  break;
    case mode::validate: validate_main(vopt);  break;
    case mode::help: std::cout << make_man_page(cli, "mantis"); break;
    }
  } else {
    auto b = res.begin();
    auto e = res.end();
    if (std::distance(b,e) > 0) {
      if (b->arg() == "build") {
        std::cout << make_man_page(build_mode, "mantis");
      } else if (b->arg() == "query") {
        std::cout << make_man_page(query_mode, "mantis");
      } else if (b->arg() == "dbgccmst") {
        std::cout << make_man_page(dbgccmst_mode, "mantis");
      } else if (b->arg() == "validate") {
        std::cout << make_man_page(validate_mode, "mantis");
      } else {
        std::cout << "There is no command \"" << b->arg() << "\"\n";
        std::cout << usage_lines(cli, "mantis") << '\n';
      }
    } else {
      std::cout << usage_lines(cli, "mantis") << '\n';
    }
  }
  
  return 0;
}
