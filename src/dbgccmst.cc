/*
 * ============================================================================
 *
 *       Filename:  dbgccmst.cc
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  2017-10-27 12:56:50 AM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Rob Johnson <robj@vmware.com>
 *   Organization:  VMware, Inc.
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

#include "ProgOpts.h"
#include "spdlog/spdlog.h"
#include "kmer.h"
#include "coloreddbg.h"
#include "common_types.h"
#include "CLI/CLI.hpp"
#include "CLI/Timer.hpp"

static void output_results(mantis::QuerySets& multi_kmers,
													 ColoredDbg<SampleObject<CQF<KeyObject>*>, KeyObject>& cdbg,
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
    }
  }
}


static void output_results_json(mantis::QuerySets& multi_kmers,
																ColoredDbg<SampleObject<CQF<KeyObject>*>, KeyObject>& cdbg,
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
      opfile << "{ \"qnum\": "       << cnt++
						 << ",  \"num_kmers\": " << kmers.size()
						 << ", \"res\": {\n";
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


int build_dbgccmst(QueryOpts& opt)
{
  std::string prefix = opt.prefix;
  std::string query_file = opt.query_file;
  std::string output_file = opt.output;
  bool use_json = opt.use_json;

  // Make sure the prefix is a full folder
  if (prefix.back() != '/') {
    prefix.push_back('/');
  }

  spdlog::logger* console = opt.console.get();
	console->info("Reading colored dbg from disk.");

	std::string cqf_file(prefix + CQF_FILE);
	std::string eqclass_file(prefix + EQCLASS_FILE);
	std::string sample_file(prefix + SAMPLEID_FILE);
	ColoredDbg<SampleObject<CQF<KeyObject>*>, KeyObject> cdbg(cqf_file,
																														eqclass_file,
																														sample_file);
  console->info("Read colored dbg with {} k-mers and {} color classes",
                cdbg.get_cqf()->size(), cdbg.get_bitvector().bit_size() / cdbg.get_num_samples());

	console->info("Reading query kmers from disk.");
	uint32_t seed = 2038074743;
  uint64_t total_kmers = 0;
  mantis::QuerySets multi_kmers = Kmer::parse_kmers(query_file.c_str(),
																										seed,
																										cdbg.range(),
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

	return EXIT_SUCCESS;
}
