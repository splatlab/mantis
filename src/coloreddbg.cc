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

#define MAX_THREADS 50
#define MAX_NUM_SAMPLES 2600
#define SAMPLE_SIZE (1ULL << 26)

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

	std::ifstream infile(opt.inlist);
	std::ifstream cutofffile(opt.cutoffs);
	//struct timeval start1, end1;
	//struct timezone tzp;
	SampleObject<CQF<KeyObject>*> *inobjects;
	CQF<KeyObject> *cqfs;


	// Allocate QF structs for input CQFs
	inobjects = (SampleObject<CQF<KeyObject>*>*)calloc(MAX_NUM_SAMPLES,
																	sizeof(SampleObject<CQF<KeyObject>*>));
	cqfs = (CQF<KeyObject>*)calloc(MAX_NUM_SAMPLES, sizeof(CQF<KeyObject>));

	// Read cutoffs files
	std::unordered_map<std::string, uint64_t> cutoffs;
	std::string sample_id;
	uint32_t cutoff;
	while (cutofffile >> sample_id >> cutoff) {
		std::pair<std::string, uint32_t> pair(last_part(sample_id, '/'), cutoff);
		cutoffs.insert(pair);
	}

	// mmap all the input cqfs
	std::string cqf_file;
	uint32_t nqf = 0;
	while (infile >> cqf_file) {
		cqfs[nqf] = CQF<KeyObject>(cqf_file, true);
		std::string sample_id = first_part(first_part(last_part(cqf_file, '/'),
																									'.'), '_');
		PRINT_CDBG("Reading CQF " << nqf << " Seed " << cqfs[nqf].seed());
		PRINT_CDBG("Sample id " << sample_id << " cut off " <<
							 cutoffs.find(sample_id)->second);
		cqfs[nqf].dump_metadata();
		inobjects[nqf] = SampleObject<CQF<KeyObject>*>(&cqfs[nqf], sample_id, nqf);
		nqf++;
	}

	std::string prefix(opt.out);
  if (prefix.back() != '/') {
    prefix += '/';
  }
  // make the output directory if it doesn't exist
  if (!mantis::fs::DirExists(prefix.c_str())) {
    mantis::fs::MakeDir(prefix.c_str());
  }

	ColoredDbg<SampleObject<CQF<KeyObject>*>, KeyObject> cdbg(opt.qbits,
																														inobjects[0].obj->keybits(),
																														inobjects[0].obj->seed(),
																														prefix, nqf);

	cdbg.build_sampleid_map(inobjects);

	PRINT_CDBG("Sampling eq classes based on " << SAMPLE_SIZE << " kmers.");
	// First construct the colored dbg on 1000 k-mers.
	cdbg_bv_map_t<__uint128_t, std::pair<uint64_t,uint64_t>> unsorted_map;

	unsorted_map = cdbg.construct(inobjects, cutoffs, unsorted_map, SAMPLE_SIZE);

	PRINT_CDBG("Number of eq classes found after sampling " <<
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

	PRINT_CDBG("Reinitializing colored DBG after the sampling phase.");
	cdbg.reinit(sorted_map);

	PRINT_CDBG("Constructing the colored dBG.");

	// Reconstruct the colored dbg using the new set of equivalence classes.
	cdbg.construct(inobjects, cutoffs, sorted_map, UINT64_MAX);

	PRINT_CDBG("Final colored dBG has " << cdbg.get_cqf()->size() <<
		" k-mers and " << cdbg.get_num_eqclasses() << " equivalence classes.");

	//cdbg.get_cqf()->dump_metadata();
	//DEBUG_CDBG(cdbg.get_cqf()->set_size());

	PRINT_CDBG("Serializing CQF and eq classes in " << prefix);
	cdbg.serialize();
	PRINT_CDBG("Serialization done.");

	return EXIT_SUCCESS;
}				/* ----------  end of function main  ---------- */
