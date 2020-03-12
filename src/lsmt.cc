/*
 * ============================================================================
 *
 *         Author:  Jamshed Khan (jamshed@cs.umd.edu)
 *   Organization:  University of Maryland, College Park
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

#include "sparsepp/spp.h"
#include "tsl/sparse_map.h"

#include "MantisFS.h"
#include "ProgOpts.h"
#include "coloreddbg.h"
#include "cdBG_merger.h"
#include "squeakrconfig.h"
#include "mantisconfig.hpp"
#include "lsmt.h"
#include "kmer.h"



int lsmt_initialize_main(LSMT_InitializeOpts &opt)
{
	spdlog::logger *console = opt.console.get();

	std::string dir = opt.dir;
	if(dir.back() != '/')	// Make sure it is a full directory.
		dir += '/';

	// Make the LSM Tree directory if it doesn't exist.
	if(!mantis::fs::DirExists(dir.c_str()))
		mantis::fs::MakeDir(dir.c_str());
	
	// Check to see if the directory exists now.
	if(!mantis::fs::DirExists(dir.c_str()))
	{
		console -> error("LSM Tree directory {} could not be created.", dir);
		exit(1);
	}

	std::ofstream pendingSamples(dir + mantis::PENDING_SAMPLES_LIST);
	pendingSamples.close();


	// Record the LSM tree parameters.
	nlohmann::json paramInfo;
	{
		std::ofstream jfile(dir + mantis::PARAM_FILE);
		if (jfile.is_open())
		{
			paramInfo["dir"] = dir;
			paramInfo["scalingFactor"] = opt.scalingFactor;
			paramInfo["kmerThreshold"] = opt.kmerThreshold;
			paramInfo["sampleThreshold"] = opt.sampleThreshold;
			paramInfo["levels"] = 0;
			paramInfo["sampleCount"] = 0;
			paramInfo["qBitInitBuild"] = opt.qBitInitBuild;


			jfile << paramInfo.dump(4);
		}
		else
		{
			console -> error("Could not write to output directory {}", dir);
			exit(1);
		}

		console -> info("LSM-tree parameters recorded at file {}.", dir + mantis::PARAM_FILE);
		
		jfile.close();
	}

	return EXIT_SUCCESS;
}



int lsmt_update_main(LSMT_UpdateOpts &opt)
{
	using LSMT_t = LSMT<SampleObject<CQF<KeyObject> *>, KeyObject>;

	spdlog::logger *console = opt.console.get();


	std::string dir = opt.dir;
	if(dir.back() != '/')	// Make sure it is a full directory.
		dir += '/';

	// Make sure the input directory exists.
	if(!mantis::fs::DirExists(dir.c_str()))
	{
		console -> error("LSM tree directory {} does not exist.", dir);
		exit(1);
	}

	
	if(!LSMT_t::is_valid_LSMT(dir, console))
		exit(1);

	std::ifstream inputList(opt.inputList);
	std::vector<std::string> inputSamples;

	if(!inputList.is_open())
	{
		console -> error("Input file {} does not exist or could not be opened.", opt.inputList);
		exit(1);
	}

	std::string sample;
	while(inputList >> sample)
		inputSamples.push_back(sample);
	
	inputList.close();

	
	LSMT_t lsmt(dir);
	lsmt.set_console(opt.console);
	lsmt.print_config();

	lsmt.update(inputSamples, opt.threadCount);


	return EXIT_SUCCESS;
}



int lsmt_query_main(LSMT_QueryOpts &opt)
{
	using LSMT_t = LSMT<SampleObject<CQF<KeyObject> *>, KeyObject>;

	spdlog::logger *console = opt.console.get();


	std::string dir = opt.dir;
	if(dir.back() != '/')	// Make sure it is a full directory.
		dir += '/';

	// Make sure the input directory exists.
	if(!mantis::fs::DirExists(dir.c_str()))
	{
		console -> error("LSM tree directory {} does not exist.", dir);
		exit(1);
	}

	
	if(!LSMT_t::is_valid_LSMT(dir, console))
		exit(1);


	console -> info("Loading LSM-tree metadata.");
	
	LSMT_t lsmt(dir);
	lsmt.set_console(opt.console);
	lsmt.print_config();


	std::string queryFile = opt.queryFile;
	std::string output = opt.output;

	uint64_t kmerLen = lsmt.kmer_len();
	console -> info("Loading query k-mers ({}-mers) from disk.", kmerLen);

	uint32_t seed = 2038074743;
	uint64_t totalKmers = 0;
	std::unordered_map<uint64_t, uint64_t> uniqueKmers;
	std::vector<std::unordered_set<uint64_t>> kmerSets = Kmer::parse_kmers(queryFile.c_str(), kmerLen,
																			totalKmers, opt.process_in_bulk,
																			uniqueKmers);
	console -> info("Total k-mers to query: {}", totalKmers);

	lsmt.query(kmerSets, output);

	return EXIT_SUCCESS;
}