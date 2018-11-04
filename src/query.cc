/*
 * ============================================================================
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
#include "mantisconfig.hpp"

void output_results(mantis::QuerySets& multi_kmers,
										ColoredDbg<SampleObject<CQF<KeyObject>*>, KeyObject>&
										cdbg, std::ofstream& opfile, bool is_bulk,
                    std::unordered_map<mantis::KmerHash, uint64_t> &uniqueKmers) {
  mantis::QueryResults qres;
	uint32_t cnt= 0;
  {
    CLI::AutoTimer timer{"Query time ", CLI::Timer::Big};
    if (is_bulk) {
        std::unordered_map<uint64_t, std::vector<uint64_t>> result = cdbg.find_samples(uniqueKmers);
        for (auto& kmers : multi_kmers) {
            opfile <<  cnt++ << '\t' << kmers.size() << '\n';
            std::vector<uint64_t> kmerCnt(cdbg.get_num_samples());
            for (auto &k : kmers) {
                if (uniqueKmers[k]) {
                    for (auto experimentId : result[k]) {
                        kmerCnt[experimentId]++;
                    }
                }
            }
            for (auto i=0; i<kmerCnt.size(); ++i) {
                if (kmerCnt[i] > 0)
                    opfile << cdbg.get_sample(i) << '\t' << kmerCnt[i] << '\n';
            }
            //++qctr;
        }
    } else {
        for (auto &kmers : multi_kmers) {
            //std::sort(kmers.begin(), kmers.end());
            opfile << cnt++ << '\t' << kmers.size() << '\n';
            mantis::QueryResult result = cdbg.find_samples(kmers);
            for (auto i = 0; i < result.size(); ++i) {
                if (result[i] > 0)
                    opfile << cdbg.get_sample(i) << '\t' << result[i] << '\n';
            }
            //++qctr;
        }
    }
  }
}

void output_results_json(mantis::QuerySets& multi_kmers,
												 ColoredDbg<SampleObject<CQF<KeyObject>*>, KeyObject>&
												 cdbg, std::ofstream& opfile, bool is_bulk,
                         std::unordered_map<mantis::KmerHash, uint64_t> &uniqueKmers) {
  mantis::QueryResults qres;
	uint32_t cnt= 0;
  {
    CLI::AutoTimer timer{"Query time ", CLI::Timer::Big};

    opfile << "[\n";
      size_t qctr{0};
      size_t nquery{multi_kmers.size()};

      if (is_bulk) {
          std::unordered_map<uint64_t, std::vector<uint64_t>> result = cdbg.find_samples(uniqueKmers);
          for (auto& kmers : multi_kmers) {
              opfile << "{ \"qnum\": " << cnt++ << ",  \"num_kmers\": " << kmers.size() << ", \"res\": {\n";
              std::vector<uint64_t> kmerCnt(cdbg.get_num_samples());
              for (auto &k : kmers) {
                  if (uniqueKmers[k]) {
                      for (auto experimentId : result[k]) {
                          kmerCnt[experimentId]++;
                      }
                  }
              }
              for (auto i=0; i<kmerCnt.size(); ++i) {
                  if (kmerCnt[i] > 0) {
                          opfile << " \"" << cdbg.get_sample(i) << "\": " << kmerCnt[i];
                      if (i != kmerCnt.size()) {
                          opfile << ",\n";
                      }
                  }
              }
              opfile << "}}";
              if (qctr < nquery - 1) { opfile << ","; }
              opfile << "\n";
              ++qctr;
          }
      }
      else {
          for (auto &kmers : multi_kmers) {
              //std::sort(kmers.begin(), kmers.end());
              opfile << "{ \"qnum\": " << cnt++ << ",  \"num_kmers\": " << kmers.size() << ", \"res\": {\n";
              mantis::QueryResult result = cdbg.find_samples(kmers);
              uint64_t sampleCntr = 0;
              for (auto it = result.begin(); it != result.end(); ++it) {
                  if (*it > 0)
                      opfile << " \"" << cdbg.get_sample(sampleCntr) << "\": " << *it;
                  if (std::next(it) != result.end()) {
                      opfile << ",\n";
                  }
                  sampleCntr++;
              }
              opfile << "}}";
              if (qctr < nquery - 1) { opfile << ","; }
              opfile << "\n";
              ++qctr;
          }
      }
    opfile << "]\n";
  }
}


/* 
 * ===  FUNCTION  =============================================================
 *         Name:  main
 *  Description:  
 * ============================================================================
 */
int query_main (QueryOpts& opt)
{
  //CLI::App app("Mantis query");

  std::string prefix = opt.prefix;
  std::string query_file = opt.query_file;
  std::string output_file = opt.output;//{"samples.output"};
  bool use_json = opt.use_json;

  // Make sure the prefix is a full folder
  if (prefix.back() != '/') {
    prefix.push_back('/');
  }
  // make the output directory if it doesn't exist
  if (!mantis::fs::DirExists(prefix.c_str())) {
    mantis::fs::MakeDir(prefix.c_str());
  }

  spdlog::logger* console = opt.console.get();
	console->info("Reading colored dbg from disk.");

	std::string dbg_file(prefix + mantis::CQF_FILE);
	std::string sample_file(prefix + mantis::SAMPLEID_FILE);
	std::vector<std::string> eqclass_files = mantis::fs::GetFilesExt(prefix.c_str(),
                                                                   mantis::EQCLASS_FILE);

	ColoredDbg<SampleObject<CQF<KeyObject>*>, KeyObject> cdbg(dbg_file,
																														eqclass_files,
																														sample_file,
																														MANTIS_DBG_IN_MEMORY);
	uint64_t kmer_size = cdbg.get_cqf()->keybits() / 2;
  console->info("Read colored dbg with {} k-mers and {} color classes",
                cdbg.get_cqf()->dist_elts(), cdbg.get_num_bitvectors());

	//cdbg.get_cqf()->dump_metadata(); 
	//CQF<KeyObject> cqf(query_file, false);
	//CQF<KeyObject>::Iterator it = cqf.begin(1);
	//mantis::QuerySet input_kmers;
	//do {
		//KeyObject k = *it;
		//input_kmers.insert(k.key);
		//++it;
	//} while (!it.done());

	//mantis::QuerySets multi_kmers;
	//multi_kmers.push_back(input_kmers);

	console->info("Reading query kmers from disk.");
	uint32_t seed = 2038074743;
	uint64_t total_kmers = 0;
    std::unordered_map<mantis::KmerHash, uint64_t> uniqueKmers;
	mantis::QuerySets multi_kmers = Kmer::parse_kmers(query_file.c_str(),
																										kmer_size,
																										total_kmers,
																										opt.process_in_bulk,
																										uniqueKmers);
	console->info("Total k-mers to query: {}", total_kmers);

	std::ofstream opfile(output_file);
	console->info("Querying the colored dbg.");

  if (use_json) {
    output_results_json(multi_kmers, cdbg, opfile, opt.process_in_bulk, uniqueKmers);
  } else {
    output_results(multi_kmers, cdbg, opfile, opt.process_in_bulk, uniqueKmers);
  }
	//std::cout << "Writing samples and abundances out." << std::endl;
	opfile.close();
	console->info("Writing done.");

	return EXIT_SUCCESS;
}				/* ----------  end of function main  ---------- */
