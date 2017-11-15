/*
 * ============================================================================
 *       Filename:  mantis.cc
 * ============================================================================
 */

#include <iostream>
#include <string>
#include <vector>

#include "ProgOpts.h"
#include "clipp.h"

#define MAX_NUM_SAMPLES 2600
#define SAMPLE_SIZE (1ULL << 26)

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

  BuildOpts bopt;
  QueryOpts qopt;
  ValidateOpts vopt;

  auto build_mode = (
                     command("build").set(selected, mode::build),
                     required("--input-list", "-i") & value("input_list", bopt.inlist) % "file containing list of input filters",
                     required("--cutoff-list", "-c") & value("cutoff_list", bopt.cutoffs) % "file containing list of experiment-specific cutoffs",
                     required("--output", "-o") & value("build_output", bopt.out) % "directory where results should be written"
                     );
  auto query_mode = (
                     command("query").set(selected, mode::query),
                     required("--input-prefix", "-i") & value("query_prefix", qopt.prefix) % "Prefix of input files.",
                     option("--output", "-o") & value("output_file", qopt.output) % "Where to write query output.",
                     option("--json","-j").set(qopt.use_json) % "Write the output in JSON format",
                     value("query", qopt.query_file) % "Prefix of input files."
                     );

  auto validate_mode = (
                     command("validate").set(selected, mode::build),
                     required("--input-list", "-i") & value("input_list", vopt.inlist) % "file containing list of input filters",
                     required("--cutoff-list", "-c") & value("cutoff_list", vopt.cutoffs) % "file containing list of experiment-specific cutoffs",
                     required("--input-prefix", "-i") & value("dbg_prefix", vopt.prefix) % "Directory containing the mantis dbg.",
                     value("query", vopt.query_file) % "Query file."
                     );

  auto cli = (
              (build_mode | query_mode | validate_mode | command("help").set(selected,mode::help) ),
              option("-v", "--version").call([]{std::cout << "version 1.0\n\n";}).doc("show version")  );

  auto res = parse(argc, argv, cli);

  //explore_options_verbose(res);

  if(res) {
    switch(selected) {
    case mode::build: build_main(bopt); /* ... */ break;
    case mode::query: query_main(qopt); /* ... */ break;
    case mode::validate: validate_main(vopt); /* ... */ break;
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
