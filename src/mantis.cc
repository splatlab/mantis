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
//#include "CLI/CLI.hpp"

#define MAX_NUM_SAMPLES 2600
#define SAMPLE_SIZE (1ULL << 26)

uint64_t start_time;

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

int query_main (QueryOpts& opt);
int build_main (BuildOpts& opt);
int validate_main (ValidateOpts& opt);

/*
 * ===  FUNCTION  =============================================================
 *         Name:  main
 *  Description:
 * ============================================================================
 */
int main ( int argc, char *argv[] ) {
  using namespace clipp;
  enum class mode {build, query, validate, help};
  mode selected = mode::help;

  auto console = spdlog::stdout_color_mt("mantis_console");

  BuildOpts bopt;
  QueryOpts qopt;
  ValidateOpts vopt;
  bopt.console = console;
  qopt.console = console;
  vopt.console = console;

	start_time = time(NULL);

  /*
  bool print_version{false};
  CLI::App app{"Mantis"};
  app.add_flag("-v,--version", print_version, "print version info");

  auto build_app = app.add_subcommand("build", "build the mantis index");
  auto query_app = app.add_subcommand("query", "query the mantis index");
  auto validate_app = app.add_subcommand("validate", "validate the mantis index");

  build_app->add_option("-i,--input-list", bopt.inlist, "file containing list of input filters")->required()->check(CLI::ExistingFile);
  build_app->add_option("-c,--cutoff-list", bopt.cutoffs, "file containing list of experiment-specific cutoffs")->required()->check(CLI::ExistingFile);
  build_app->add_option("-o,--output", bopt.out, "directory where results should be written")->required();

  query_app->add_flag("-j,--json", qopt.use_json, "Write the output in JSON format");
  query_app->add_option("-p,--input-prefix", qopt.prefix, "Prefix of input files.")->required()->check(CLI::ExistingDirectory);
  query_app->add_option("-o,--output", qopt.output, "Where to write query output.");
  query_app->add_option("query", qopt.query_file, "Prefix of input files.")->required()->check(CLI::ExistingFile);

  validate_app->add_option("-i,--input-list", vopt.inlist, "file containing list of input filters")->required()->check(CLI::ExistingFile);
  validate_app->add_option("-c,--cutoff-list", vopt.cutoffs, "file containing list of experiment-specific cutoffs")->required()->check(CLI::ExistingFile);
  validate_app->add_option("-p,--input-prefix", vopt.prefix, "Directory containing the mantis dbg.")->required()->check(CLI::ExistingDirectory);
  validate_app->add_option("query", vopt.query_file,"Query file.")->required()->check(CLI::ExistingFile);

  CLI11_PARSE(app, argc, argv);

  if (print_version) {
    std::cerr << "mantis " << 0.1 << "\n";
    return 0;
  }

  if (app.got_subcommand(build_app)) {
    return build_main(bopt);
  } else if (app.got_subcommand(query_app)) {
    return query_main(qopt);
  } else if (app.got_subcommand(validate_app)) {
    return validate_main(vopt);
  } else {
    std::cout << "I don't know the requested sub-command\n";
    return 1;
  }
  */

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
                     required("-i", "--input-list") & value(ensure_file_exists, "input_list", bopt.inlist) % "file containing list of input filters",
                     required("-o", "--output") & value("build_output", bopt.out) % "directory where results should be written"
                     );
  auto query_mode = (
                     command("query").set(selected, mode::query),
                     option("-j", "--json").set(qopt.use_json) % "Write the output in JSON format",
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

  auto cli = (
              (build_mode | query_mode | validate_mode | command("help").set(selected,mode::help) ),
              option("-v", "--version").call([]{std::cout << "version 1.0\n\n";}).doc("show version")  );

  assert(build_mode.flags_are_prefix_free());
  assert(query_mode.flags_are_prefix_free());
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
