/*
 * =====================================================================================
 *
 *         Author:  Prashant Pandey (), ppandey@cs.stonybrook.edu
 *   Organization:  Stony Brook University
 *
 * =====================================================================================
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

#include "MantisFS.h"
#include "ProgOpts.h"
#include "kmer.h"
#include "coloreddbg.h"
#include "mantisconfig.hpp"

#include	<stdlib.h>

/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  main
 *  Description:  
 * =====================================================================================
 */
	int
validate_main ( ValidateOpts& opt )
{

	spdlog::logger* console = opt.console.get();
	// Read experiment CQFs
	std::ifstream infile(opt.inlist);
  uint64_t num_samples{0};
  if (infile.is_open()) {
    std::string line;
    while (std::getline(infile, line)) { ++num_samples; }
    infile.clear();
    infile.seekg(0, std::ios::beg);
  } else {
    console->error("Input filter list {} does not exist or could not be opened.", opt.inlist);
    std::exit(1);
  }
	SampleObject<CQF<KeyObject>*> *inobjects;
	CQF<KeyObject> *cqfs;

	// Allocate QF structs for input CQFs
	inobjects = (SampleObject<CQF<KeyObject>*>*)calloc(num_samples,
																										 sizeof(SampleObject<CQF<KeyObject>*>));
	cqfs = (CQF<KeyObject>*)calloc(num_samples, sizeof(CQF<KeyObject>));


	// mmap all the input cqfs
	std::string cqf_file;
	uint32_t nqf = 0;
	while (infile >> cqf_file) {
		if (!mantis::fs::FileExists(cqf_file.c_str())) {
			console->error("Squeakr file {} does not exist.", cqf_file);
			exit(1);
		}
		cqfs[nqf] = CQF<KeyObject>(cqf_file, CQF_FREAD);
		std::string sample_id = first_part(first_part(last_part(cqf_file, '/'),
																									'.'), '_');
		console->info("Reading CQF {} Seed {}", nqf, cqfs[nqf].seed());
		console->info("Sample id {} cut-off {}", sample_id);
		cqfs[nqf].dump_metadata();
		inobjects[nqf] = SampleObject<CQF<KeyObject>*>(&cqfs[nqf],
																									 sample_id, nqf);
		nqf++;
	}

	std::string prefix = opt.prefix;
	if (prefix.back() != '/') {
		prefix += '/';
	}
	// make the output directory if it doesn't exist
	if (!mantis::fs::DirExists(prefix.c_str())) {
		mantis::fs::MakeDir(prefix.c_str());
	}

	// Read the colored dBG
	console->info("Reading colored dbg from disk.");
	std::string dbg_file(prefix + mantis::CQF_FILE);
	std::string sample_file(prefix + mantis::SAMPLEID_FILE);
	std::vector<std::string> eqclass_files = mantis::fs::GetFilesExt(prefix.c_str(),
																																	 mantis::EQCLASS_FILE);

	ColoredDbg<SampleObject<CQF<KeyObject>*>, KeyObject> cdbg(dbg_file,
																														eqclass_files,
																														sample_file);

	uint64_t kmer_size = cdbg.get_cqf()->keybits() / 2;
	console->info("Read colored dbg with {} k-mers and {} color classes",
								cdbg.get_cqf()->dist_elts(), cdbg.get_num_bitvectors());

	std::string query_file = opt.query_file;
	console->info("Reading query kmers from disk.");
	uint32_t seed = 2038074743;
	uint64_t total_kmers = 0;
	mantis::QuerySets multi_kmers = Kmer::parse_kmers(query_file.c_str(),
																										kmer_size,
																										total_kmers);
	console->info("Total k-mers to query: {}", total_kmers);

	// Query kmers in each experiment CQF
	// Maintain the fraction of kmers present in each experiment CQF.
	std::vector<std::unordered_map<uint64_t, float>> ground_truth;
	std::vector<std::unordered_map<uint64_t, uint64_t>> cdbg_output;
	bool fail{false};
	for (auto kmers : multi_kmers) {
		std::unordered_map<uint64_t, float> fraction_present;
		for (uint64_t i = 0; i < nqf; i++) {
			for (auto kmer : kmers) {
				KeyObject k(kmer, 0, 0);
				uint64_t count = cqfs[i].query(k, 0);
				fraction_present[inobjects[i].id] += 1;
			}
		}
		// Query kmers in the cdbg
		std::unordered_map<uint64_t, uint64_t> result = cdbg.find_samples(kmers);

		// Validate the cdbg output
		for (uint64_t i = 0; i < nqf; i++)
			if (fraction_present[i] != result[i]) {
				console->info("Failed for sample: {} original CQF {} cdbg {}",
											inobjects[i].sample_id, fraction_present[i], result[i]);
				fail = true;
				//abort();
			}
		ground_truth.push_back(fraction_present);
		cdbg_output.push_back(result);
	}
	if (fail)
		console->info("Mantis validation failed!");
	else
		console->info("Mantis validation passed!");

#if 0
	// This is x-axis
	// For one query set
	// Divide aggregate by the total kmers in the query to get the fraction.
	uint64_t cnt = 0;
	std::unordered_map<uint64_t, float> fraction_present = ground_truth[0];
	for (uint64_t i = 0; i < nqf; i++)
		fraction_present[i] = (fraction_present[i] / (float)multi_kmers[cnt++].size()) * 100;

	std::vector<std::vector<std::string>> buckets;
	for (auto it : fraction_present)
		buckets[it.second / 5].push_back(cdbg.get_sample(it.first));
#endif

	return EXIT_SUCCESS;
}				/* ----------  end of function main  ---------- */

