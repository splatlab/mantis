/*
 * ============================================================================
 *
 *       Filename:  main.cc
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  2016-11-10 03:31:54 PM
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
#include "MantisFS.h"
#include "sparsepp/spp.h"
#include "tsl/sparse_map.h"
#include "ProgOpts.h"
#include "coloreddbg.h"
#include "json.hpp"
#include "mantis_utils.hpp"
#include "mantisconfig.hpp"

// This function read one byte from each page in the iterator buffer.
uint64_t tmp_sum;
void handler_function(union sigval sv) {
	CQF<KeyObject>::Iterator& it(*((CQF<KeyObject>::Iterator*)sv.sival_ptr));
	unsigned char *start = (unsigned char*)(it.iter.qf->metadata) + it.last_prefetch_offset;
	unsigned char *counter = (unsigned char*)(it.iter.qf->metadata) + it.last_prefetch_offset;
	for (;counter < start + it.buffer_size; counter += 4096) {
		tmp_sum += *counter;
	}
}

/*
 * ===  FUNCTION  =============================================================
 *         Name:  main
 *  Description:  
 * ============================================================================
 */
	int
build_main ( BuildOpts& opt )
{
	/* calling asyc read init */
	struct aioinit aioinit;
	memset(&aioinit, 0, sizeof(struct aioinit));
	aioinit.aio_num = 2500;
	aioinit.aio_threads = 100;
	aioinit.aio_idle_time = 60;
	aio_init(&aioinit);

	spdlog::logger* console = opt.console.get();
	std::ifstream infile(opt.inlist);
  uint64_t num_samples{0};
  if (infile.is_open()) {
    std::string line;
    while (std::getline(infile, line)) { ++num_samples; }
    infile.clear();
    infile.seekg(0, std::ios::beg);
    console->info("Will build mantis index over {} input experiments.", num_samples);
  } else {
    console->error("Input file {} does not exist or could not be opened.", opt.inlist);
    std::exit(1);
  }
	//struct timeval start1, end1;
	//struct timezone tzp;

  // C++-izing
  // This is C++ ... not C.  Why do we have raw pointers
  // and calloc them.  We use vectors instead.
	//SampleObject<CQF<KeyObject>*> *inobjects;
	//CQF<KeyObject> *cqfs;
  std::vector<SampleObject<CQF<KeyObject>*>> inobjects;
  std::vector<CQF<KeyObject>> cqfs;

  /** try and create the output directory
   *  and write a file to it.  Complain to the user
   *  and exit if we cannot.
   **/
	std::string prefix(opt.out);
	if (prefix.back() != '/') {
		prefix += '/';
	}
	// make the output directory if it doesn't exist
	if (!mantis::fs::DirExists(prefix.c_str())) {
		mantis::fs::MakeDir(prefix.c_str());
	}
	// check to see if the output dir exists now
	if (!mantis::fs::DirExists(prefix.c_str())) {
		console->error("Output dir {} could not be successfully created.", prefix);
		exit(1);
	}

  // If we made it this far, record relevant meta information in the output directory
  nlohmann::json minfo;
  {
    std::ofstream jfile(prefix + "/" + mantis::meta_file_name);
    if (jfile.is_open()) {
      minfo = opt.to_json();
      minfo["start_time"] = mantis::get_current_time_as_string();
      minfo["mantis_version"] = mantis::version;
      minfo["index_version"] = mantis::index_version;
      jfile << minfo.dump(4);
    } else {
      console->error("Could not write to output directory {}", prefix);
      exit(1);
    }
    jfile.close();
  }

	// reserve QF structs for input CQFs
  inobjects.reserve(num_samples);
  cqfs.reserve(num_samples);

	// mmap all the input cqfs
	std::string cqf_file;
	uint32_t nqf = 0;
	uint32_t cutoff;
	console->info("Reading input Squeakr files.");
	while (infile >> cqf_file >> cutoff) {
		if (!mantis::fs::FileExists(cqf_file.c_str())) {
			console->error("Squeakr file {} does not exist.", cqf_file);
			exit(1);
		}
    cqfs.emplace_back(cqf_file, true);
		std::string sample_id = first_part(first_part(last_part(cqf_file, '/'),
																									'.'), '_');
		console->info("Reading CQF {} Seed {}",nqf, cqfs[nqf].seed());
		console->info("Sample id {} cut off {}", sample_id, cutoff);
		cqfs.back().dump_metadata();
    inobjects.emplace_back(&cqfs[nqf], cutoff, sample_id, nqf);
    nqf++;
	}

	ColoredDbg<SampleObject<CQF<KeyObject>*>, KeyObject> cdbg(opt.qbits,
																														inobjects[0].obj->keybits(),
																														inobjects[0].obj->seed(),
																														prefix, nqf);
	cdbg.set_console(console);
	if (opt.flush_eqclass_dist) {
		cdbg.set_flush_eqclass_dist();
  }

	cdbg.build_sampleid_map(inobjects.data());

	console->info("Sampling eq classes based on {} kmers", mantis::SAMPLE_SIZE);
	// First construct the colored dbg on initial SAMPLE_SIZE k-mers.
	default_cdbg_bv_map_t unsorted_map;

	unsorted_map = cdbg.construct(inobjects.data(), mantis::SAMPLE_SIZE);

	console->info("Number of eq classes found after sampling {}",
								unsorted_map.size());

	// Sort equivalence classes based on their abundances.
	std::multimap<uint64_t, __uint128_t, std::greater<uint64_t>> sorted;
	for (auto& it : unsorted_map) {
		//DEBUG_CDBG(it.second.second << " " << it.second.first << " " << it.first.data());
		sorted.insert(std::pair<uint64_t, __uint128_t>(it.second.second,
																									 it.first));
		//sorted[it.second.second] = it.first;
	}
	cdbg_bv_map_t<__uint128_t, std::pair<uint64_t,uint64_t>> sorted_map;
	//DEBUG_CDBG("After sorting.");
	uint64_t i = 1;
	for (auto& it : sorted) {
		//DEBUG_CDBG(it.first << " " << it.second.data());
		std::pair<uint64_t, uint64_t> val(i, 0);
		std::pair<__uint128_t, std::pair<uint64_t, uint64_t>> keyval(it.second, val);
		sorted_map.insert(keyval);
		i++;
	}

	console->info("Reinitializing colored DBG after the sampling phase.");
	cdbg.reinit(sorted_map);

	console->info("Constructing the colored dBG.");

	// Reconstruct the colored dbg using the new set of equivalence classes.
	cdbg.construct(inobjects.data(), std::numeric_limits<uint64_t>::max());

	console->info("Final colored dBG has {} k-mers and {} equivalence classes",
								cdbg.get_cqf()->size(), cdbg.get_num_eqclasses());

	//cdbg.get_cqf()->dump_metadata();
	//DEBUG_CDBG(cdbg.get_cqf()->set_size());

	console->info("Serializing CQF and eq classes in {}", prefix);
	cdbg.serialize();
	console->info("Serialization done.");

  {
    std::ofstream jfile(prefix + "/" + mantis::meta_file_name);
    if (jfile.is_open()) {
      minfo["end_time"] = mantis::get_current_time_as_string();
      jfile << minfo.dump(4);
    } else {
      console->error("Could not write to output directory {}", prefix);
    }
    jfile.close();
  }

  return EXIT_SUCCESS;
}				/* ----------  end of function main  ---------- */
