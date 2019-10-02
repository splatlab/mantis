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
#include "CdBG_Merger.h"
#include "squeakrconfig.h"
#include "mantisconfig.hpp"
#include "lsmt.h"



void lsmt_initialize_main(LSMT_InitializeOpts &opt)
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
		console -> error("LSM Tree dir {} could not be created.", dir);
		exit(1);
	}

	const char PARAM_FILE[] = "lsmt-params.json"; // TODO: move it


	// Record the LSM tree parameters.
	nlohmann::json paramInfo;
	{
		std::ofstream jfile(dir + PARAM_FILE);
		if (jfile.is_open())
		{
			paramInfo["dir"] = dir;
			paramInfo["scalingFactor"] = opt.scalingFactor;
			paramInfo["kmerThreshold"] = opt.kmerThreshold;
			paramInfo["sampleThreshold"] = opt.sampleThreshold;
			paramInfo["levels"] = 0;


			jfile << paramInfo.dump(4);
		}
		else
		{
			console -> error("Could not write to output directory {}", dir);
			exit(1);
		}

		console -> info("LSM-tree parameters recorded at file {}.", dir + PARAM_FILE);
		
		jfile.close();
	}
}
