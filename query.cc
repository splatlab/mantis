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

#include "kmer.h"
#include "coloreddbg.h"
#include "CLI/CLI.hpp"
#include "CLI/Timer.hpp"

#define MAX_NUM_SAMPLES 2600
#define OUTPUT_FILE "samples.output"

/* 
 * ===  FUNCTION  =============================================================
 *         Name:  main
 *  Description:  
 * ============================================================================
 */
int main ( int argc, char *argv[] )
{
  CLI::App app("Mantis query");

  std::string prefix;
  std::string query_file;
  std::string output_file{"samples.output"};
  app.add_option("-i,--input-prefix", prefix, "Prefix of input files.")->required();
  app.add_option("-o,--outout", output_file, "Where to write query output.");
  app.add_option("query", query_file, "Prefix of input files.")->required();

  CLI11_PARSE(app, argc, argv);

  // Make sure the prefix is a full folder
  if (prefix.back() != '/') {
    prefix.push_back('/');
  }

	std::cout << "Reading colored dbg from disk." << std::endl;
	std::string cqf_file(prefix + CQF_FILE);
	std::string eqclass_file(prefix + EQCLASS_FILE);
	std::string sample_file(prefix + SAMPLEID_FILE);
	ColoredDbg<SampleObject<CQF<KeyObject>*>, KeyObject> cdbg(cqf_file,
																														eqclass_file,
																														sample_file);
	std::cout << "Read colored dbg with " << cdbg.get_cqf()->size() << " k-mers and "
						 << cdbg.get_bitvector().bit_size() / cdbg.get_num_samples() <<
						 " equivalence classes." << std::endl;
	
	//cdbg.get_cqf()->dump_metadata(); 
	//std::string query_file(argv[2]);
	//CQF<KeyObject> cqf(query_file);
	//CQF<KeyObject>::Iterator it = cqf.begin(1);
	//std::vector<uint64_t> input_kmers;
	//do {
		//KeyObject k = *it;
		//input_kmers.push_back(k.key);
		//++it;
	//} while (!it.done());

	std::cout << "Reading query kmers from disk." << std::endl;
	uint32_t seed = 2038074743;
	std::vector<std::unordered_set<uint64_t>> multi_kmers = Kmer::parse_kmers(query_file.c_str(),
																																		 seed,
																																		 cdbg.range());

	std::ofstream opfile(output_file);
	std::cout << "Querying the colored dbg." << std::endl;

	uint32_t cnt= 0;
  {
    CLI::AutoTimer timer{"Query time ", CLI::Timer::Big};
    for (auto kmers : multi_kmers) {
      opfile << cnt++ << "\t" << kmers.size() << std::endl;
      std::unordered_map<uint64_t, uint64_t> result = cdbg.find_samples(kmers);
      for (auto it = result.begin(); it != result.end(); ++it)
        opfile << cdbg.get_sample(it->first) << " " << it->second << std::endl;
    }
  }
	//std::cout << "Writing samples and abundances out." << std::endl;
	opfile.close();
	std::cout << "Writing done." << std::endl;

	return EXIT_SUCCESS;
}				/* ----------  end of function main  ---------- */
