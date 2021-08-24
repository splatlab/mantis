/*
 * ============================================================================
 *
 *         Author:  Jamshed Khan (jamshed@umd.edu)
 *   Organization:  University of Maryland, College Park
 *
 * ============================================================================
 */


#include "MantisFS.h"
#include "ProgOpts.h"
#include "coloreddbg.h"
#include "lsmt.h"
#include "kmer.h"
 
#include <fstream>
#include <string>
#include <vector>



int lsmt_initialize_main(LSMT_InitializeOpts& opt)
{
	spdlog::logger* console = opt.console.get();

	std::string dir = opt.dir;
	if(dir.back() != '/')	// Make sure it is a full directory.
		dir += '/';

	// Make the LSM tree directory if it doesn't exist.
	if(!mantis::fs::DirExists(dir.c_str()))
		mantis::fs::MakeDir(dir.c_str());
	
	// Check to see if the directory exists now.
	if(!mantis::fs::DirExists(dir.c_str()))
	{
		console->error("Error creating LSM tree directory {}. Aborting.", dir);
		std::exit(EXIT_FAILURE);
	}

	// Create the pending samples list file.
	std::ofstream pendingSamples(dir + mantis::PENDING_SAMPLES_LIST);
	if(!pendingSamples.is_open())
	{
		console->error("Error creating the pending samples list file {}. Aborting.", dir + mantis::PENDING_SAMPLES_LIST);
		std::exit(EXIT_FAILURE);
	}
	pendingSamples.close();


	// Record the LSM tree parameters.
	nlohmann::json paramInfo;
	{
		std::ofstream jfile(dir + mantis::PARAM_FILE);
		if (jfile.is_open())
		{
			paramInfo["dir"] = dir;
			paramInfo["scaling_factor"] = opt.scaling_factor;
			paramInfo["kmer_threshold"] = opt.kmer_threshold;
			paramInfo["cqf_count_threshold"] = opt.cqf_count_threshold;
			paramInfo["sample_threshold"] = opt.sample_threshold;
			paramInfo["levels"] = 0;
			paramInfo["sample_count"] = 0;
			paramInfo["qBit_init_build"] = opt.qBit_init_build;


			jfile << paramInfo.dump(4);
		}
		else
		{
			console->error("Error opening the LSM tree parameters file {}. Aborting.", dir + mantis::PARAM_FILE);
			std::exit(EXIT_FAILURE);
		}

		console->info("LSM tree parameters recorded at file {}.", dir + mantis::PARAM_FILE);
		
		jfile.close();
	}


	return EXIT_SUCCESS;
}



int lsmt_update_main(LSMT_UpdateOpts& opt)
{
	using LSMT_t = LSMT<SampleObject<CQF<KeyObject> *>, KeyObject>;

	spdlog::logger* console = opt.console.get();


	std::string dir = opt.dir;
	if(dir.back() != '/')	// Make sure it is a full directory.
		dir += '/';

	// Make sure the input directory exists.
	if(!mantis::fs::DirExists(dir.c_str()))
	{
		console->error("LSM tree directory {} does not exist.", dir);
		std::exit(EXIT_FAILURE);
	}


	if(!LSMT_t::is_valid_LSMT(dir, console))
	{
		console->error("LSM tree structure is not valid. Aborting.");
		std::exit(EXIT_FAILURE);
	}

	std::ifstream input_list(opt.input_list);
	std::vector<std::string> input_samples;

	if(!input_list.is_open())
	{
		console->error("Input file {} does not exist or could not be opened.", opt.input_list);
		std::exit(EXIT_FAILURE);
	}

	std::string sample;
	while(input_list >> sample)
		input_samples.push_back(sample);
	
	input_list.close();

	
	LSMT_t lsmt(dir);
	lsmt.set_console(opt.console);
	lsmt.print_config();

	lsmt.update(input_samples, opt.thread_count);


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


int lsmt_query_main(QueryOpts& opt)
{
	typedef LSMT<SampleObject<CQF<KeyObject> *>, KeyObject> LSMT_t;

	spdlog::logger* const logger = opt.console.get();

	if(!opt.process_in_bulk)
	{
		logger->error("Only qurying in bulk is supported currently. Aborting.\n");
		std::exit(EXIT_FAILURE);
	}

	const std::string lsmt_dir(opt.prefix + "/");
	if(!mantis::fs::DirExists(lsmt_dir.c_str()))
	{
		logger->error("LSM-tree directory does not exist. Aborting.");
		std::exit(EXIT_FAILURE);
	}

	if(!LSMT_t::is_valid_LSMT(lsmt_dir, logger))
	{
		logger->error("LSM-tree is not valid. Aborting.\n");
		std::exit(EXIT_FAILURE);
	}



	LSMT_t lsmt(lsmt_dir);
	lsmt.set_console(opt.console);
	lsmt.print_config();


	lsmt.query(opt);


	return EXIT_SUCCESS;
}
