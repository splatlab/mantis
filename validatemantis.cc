/*
 * =====================================================================================
 *
 *       Filename:  validate_mantis.cc
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  2017-10-31 12:13:49 PM
 *       Revision:  none
 *       Compiler:  gcc
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

#include "kmer.h"
#include "coloreddbg.h"

#define MAX_NUM_SAMPLES 2600
#define SAMPLE_SIZE (1ULL << 26)

#include	<stdlib.h>

/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  main
 *  Description:  
 * =====================================================================================
 */
	int
main ( int argc, char *argv[] )
{
	if (argc < 5) {
		std::cout << "Not suffcient args." << std::endl;
		abort();
	}

	// Read experiment CQFs and cutoffs
	std::ifstream infile(argv[1]);
	std::ifstream cutofffile(argv[2]);
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
		cqfs[nqf] = CQF<KeyObject>(cqf_file, false);
		std::string sample_id = first_part(first_part(last_part(cqf_file, '/'),
																									'.'), '_');
		std::cout << "Reading CQF " << nqf << " Seed " << cqfs[nqf].seed() <<
			std::endl;
		std::cout << "Sample id " << sample_id << " cut off " <<
			cutoffs.find(sample_id)->second << std::endl;
		cqfs[nqf].dump_metadata();
		inobjects[nqf] = SampleObject<CQF<KeyObject>*>(&cqfs[nqf], sample_id, nqf);
		nqf++;
	}

	// Read the colored dBG
	std::string prefix(argv[3]);
	std::cout << "Reading colored dbg from disk." << std::endl;
	std::string dbg_file(prefix + CQF_FILE);
	std::string eqclass_file(prefix + EQCLASS_FILE);
	std::string sample_file(prefix + SAMPLEID_FILE);
	ColoredDbg<SampleObject<CQF<KeyObject>*>, KeyObject> cdbg(dbg_file,
																														eqclass_file,
																														sample_file);
	std::cout << "Read colored dbg with " << cdbg.get_cqf()->size() << " k-mers and "
		<< cdbg.get_bitvector().bit_size() / cdbg.get_num_samples() <<
		" equivalence classes." << std::endl;

	// Read query k-mers.
	std::cout << "Reading query kmers from disk." << std::endl;
	uint32_t seed = 2038074743;
	std::vector<std::unordered_set<uint64_t>> multi_kmers = Kmer::parse_kmers(argv[4],
																																		 seed,
																																		 cdbg.range());

	//typename CQF<KeyObject>::Iterator it = inobjects[0].obj->begin(0);
	//std::ofstream kmerlist("kmerlist.dump");
	//while(!it.done()) {
		//KeyObject k = *it;
		//kmerlist << int_to_str(HashUtil::hash_64i(k.key, BITMASK(2*K))) <<
			//"\t" << k.count << std::endl;
		//++it;
	//}
	//kmerlist.close();

	// Query kmers in each experiment CQF ignoring kmers below the cutoff.
	// Maintain the fraction of kmers present in each experiment CQF.
	std::vector<std::unordered_map<uint64_t, float>> ground_truth;
	std::vector<std::unordered_map<uint64_t, uint64_t>> cdbg_output;
	for (auto kmers : multi_kmers) {
		std::unordered_map<uint64_t, float> fraction_present;
		for (uint64_t i = 0; i < nqf; i++) {
			uint32_t cutoff;
			auto it = cutoffs.find(inobjects[i].sample_id);
			if (it == cutoffs.end()) {
				std::cerr << "Sample id " <<  inobjects[i].sample_id << " not found in"
					<< " cutoff list." << std::endl;
				abort();
			} else
				cutoff = it->second;
			for (auto kmer : kmers) {
				KeyObject k(kmer, 0, 0);
				uint64_t count = cqfs[i].query(k);
				if (count < cutoff)
					continue;
				else
					fraction_present[inobjects[i].id] += 1;
			}
		}
		// Query kmers in the cdbg
		std::unordered_map<uint64_t, uint64_t> result = cdbg.find_samples(kmers);

		// Validate the cdbg output
		for (uint64_t i = 0; i < nqf; i++)
			if (fraction_present[i] != result[i]) {
				std::cout << "Failed for sample: " << inobjects[i].sample_id << 
					" original CQF " << fraction_present[i] << " cdbg " << 
					result[i] << std::endl;
				//abort();
			}
		ground_truth.push_back(fraction_present);
		cdbg_output.push_back(result);
	}
	std::cout << "Colored dBG output is validated." << std::endl;

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

