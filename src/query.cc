/*
 * ============================================================================
 *
 *       Filename:  query.cc
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  2017-10-27 12:56:50 AM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Prashant Pandey (), ppandey@cs.stonybrook.edu
 *   Organization:  Stony Brook University
 *
 * ============================================================================
 */


#include <iostream>
#include <fstream>
#include <utility>
#include <algorithm>
#include <string>
#include <vector>
#include <map>
#include <queue>
#include <set>
#include <unordered_set>
#include <set>
#include <bitset>
#include <cassert>
#include <fstream>

#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdlib.h>
#include <fcntl.h>
#include <unistd.h>
#include <sys/resource.h>
#include <sys/stat.h>
#include <sys/time.h>
#include <sys/mman.h>
#include <openssl/rand.h>

#include "MantisFS.h"
#include "ProgOpts.h"
#include "spdlog/spdlog.h"
#include "kmer.h"
#include "coloreddbg.h"
#include "common_types.h"
#include "CLI/CLI.hpp"
#include "CLI/Timer.hpp"

#include "sockets.h"

#define MAX_NUM_SAMPLES 2600
#define OUTPUT_FILE "samples.output"

void output_results(mantis::QuerySets& multi_kmers,   ColoredDbg<SampleObject<CQF<KeyObject>*>, KeyObject>& cdbg,
                    std::ofstream& opfile) {
  mantis::QueryResults qres;
  uint32_t cnt= 0;
  {
    CLI::AutoTimer timer{"Query time ", CLI::Timer::Big};
    //size_t qctr{0};
    //size_t nquery{multi_kmers.size()};
    for (auto& kmers : multi_kmers) {
      //std::sort(kmers.begin(), kmers.end());
      opfile <<  cnt++ << '\t' << kmers.size() << '\n';
      mantis::QueryResult result = cdbg.find_samples(kmers);
      for (auto it = result.begin(); it != result.end(); ++it) {
        opfile << cdbg.get_sample(it->first) << '\t' << it->second << '\n';
      }
      //++qctr;
    }
  }
}


void output_results_json(mantis::QuerySets& multi_kmers,  ColoredDbg<SampleObject<CQF<KeyObject>*>, KeyObject>& cdbg,
                    std::ofstream& opfile) {
  mantis::QueryResults qres;
  uint32_t cnt= 0;
  {
    CLI::AutoTimer timer{"Query time ", CLI::Timer::Big};
    opfile << "[\n";
    size_t qctr{0};
    size_t nquery{multi_kmers.size()};
    for (auto& kmers : multi_kmers) {
      //std::sort(kmers.begin(), kmers.end());
      opfile << "{ \"qnum\": " << cnt++ << ",  \"num_kmers\": " << kmers.size() << ", \"res\": {\n";
      mantis::QueryResult result = cdbg.find_samples(kmers);
      for (auto it = result.begin(); it != result.end(); ++it) {
        opfile << " \"" <<cdbg.get_sample(it->first) << "\": " << it->second ;
        if (std::next(it) != result.end()) {
          opfile << ",\n";
        }
      }
      opfile << "}}";
      if (qctr < nquery - 1) { opfile << ","; }
      opfile << "\n";
      ++qctr;
    }
    opfile << "]\n";
  }

}

void run_query(std::string query_file, std::string output_file, 
               ColoredDbg<SampleObject<CQF<KeyObject>*>, KeyObject>& cdbg, 
               spdlog::logger * console, 
               bool use_json) {
  uint32_t seed = 2038074743;
  uint64_t total_kmers = 0;

  uint64_t kmer_size = cdbg.get_cqf()->keybits() / 2;
  console->info("Use colored dbg with {} k-mers and {} color classes",
                cdbg.get_cqf()->size(), cdbg.get_num_bitvectors());
  console->info("K-mer size: {}", kmer_size);
  console->info("Query file: {}", query_file);
  console->info("Output file: {}", output_file);

  mantis::QuerySets multi_kmers = Kmer::parse_kmers(query_file.c_str(),
                                                    seed,
                                                    cdbg.range(),
                                                    kmer_size,
                                                    total_kmers);
  console->info("Total k-mers to query: {}", total_kmers);

  std::ofstream opfile(output_file);
  console->info("Querying the colored dbg.");

  if (use_json) {
    output_results_json(multi_kmers, cdbg, opfile);
  } else {
    output_results(multi_kmers, cdbg, opfile);
  }

  opfile.close();
  console->info("Writing done.");
}


/* 
 * ===  FUNCTION  =============================================================
 *         Name:  query main
 *  Description:  
 * ============================================================================
 */
int query_main(QueryOpts& opt)
{
  std::string prefix = opt.prefix;
  std::string query_file = opt.query_file;
  std::string output_file = opt.output;//{"samples.output"};
  bool use_json = opt.use_json;
  // Make sure the prefix is a full folder
  if (prefix.back() != '/') {
    prefix.push_back('/');
  }

  spdlog::logger* console = opt.console.get();
  console->info("Reading colored dbg from disk.");

  std::string cqf_file(prefix + CQF_FILE);
  std::string sample_file(prefix + SAMPLEID_FILE);
  std::vector<std::string> eqclass_files = mantis::fs::GetFilesExt(prefix.c_str(),
                                                       std::string(EQCLASS_FILE).c_str());

  ColoredDbg<SampleObject<CQF<KeyObject>*>, KeyObject> cdbg(cqf_file,
                                                            eqclass_files,
                                                            sample_file);
  run_query(query_file, output_file, cdbg, console, use_json); 

  return EXIT_SUCCESS;
}       /* ----------  end of function main  ---------- */


/* 
 * ===  FUNCTION  =============================================================
 *         Name:  server main
 *  Description:  
 * ============================================================================
 */
int server_main(QueryOpts& opt)
{
  try {
    std::string prefix = opt.prefix;
    bool use_json = opt.use_json;

    // Make sure the prefix is a full folder
    if (prefix.back() != '/') {
      prefix.push_back('/');
    }

    spdlog::logger* console = opt.console.get();
    console->info("Reading colored dbg from disk: " + prefix);

    std::string cqf_file(prefix + CQF_FILE);
    std::string sample_file(prefix + SAMPLEID_FILE);
    std::vector<std::string> eqclass_files = mantis::fs::GetFilesExt(prefix.c_str(),
                                                         std::string(EQCLASS_FILE).c_str());
    ColoredDbg<SampleObject<CQF<KeyObject>*>, KeyObject> cdbg(cqf_file,
                                                              eqclass_files,
                                                              sample_file);

    boost::asio::io_service io_service;
    server s(io_service, 23901, cdbg, console, use_json);

    console->info("Reading query kmers from disk.");
    console->info("Run server accepting queries.");
    io_service.run();

    // read from stdin instead of socket
    // while (1) {
    //   std::string query_file = opt.query_file;
    //   std::string output_file = opt.output;  //{"samples.output"};

    //   cout << "Enter query file: " << endl;
    //   cin >> query_file;
    //   cout << "Enter output file: " << endl;
    //   cin >> output_file;

    //   run_query(query_file, output_file, cdbg, console, use_json); 
    // }
  }
  catch (std::exception& e) {
    std::cerr << "Exception: " << e.what() << "\n";
  }

  return EXIT_SUCCESS;
}       /* ----------  end of function main  ---------- */

