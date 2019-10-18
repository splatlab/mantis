/*
 * ============================================================================
 *
 *         Author:  Prashant Pandey (), ppandey@cs.stonybrook.edu
 * 					Jamshed Khan (), jamshed@umd.edu
 *   Organization:  Stony Brook University, University of Maryland
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
#include "squeakrconfig.h"
#include "json.hpp"
#include "mantis_utils.hpp"
#include "mantisconfig.hpp"

/*
 * ===  FUNCTION  =============================================================
 *         Name:  main
 *  Description:  
 * ============================================================================
 */
	int
build_main ( BuildOpts& opt )
{
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

	std::vector<SampleObject<CQF<KeyObject>*>> inobjects;
  std::vector<CQF<KeyObject>> cqfs;

	// reserve QF structs for input CQFs
  inobjects.reserve(num_samples);
  cqfs.reserve(num_samples);

	// mmap all the input cqfs
	std::string squeakr_file;
	uint32_t nqf = 0;
	uint32_t kmer_size{0};
	console->info("Reading input Squeakr files.");
	while (infile >> squeakr_file) {
		if (!mantis::fs::FileExists(squeakr_file.c_str())) {
			console->error("Squeakr file {} does not exist.", squeakr_file);
			exit(1);
		}
		squeakr::squeakrconfig config;
		int ret = squeakr::read_config(squeakr_file, &config);
		if (ret == squeakr::SQUEAKR_INVALID_VERSION) {
			console->error("Squeakr index version is invalid. Expected: {} Available: {}",
										 squeakr::INDEX_VERSION, config.version);
			exit(1);
		}
		if (ret == squeakr::SQUEAKR_INVALID_ENDIAN) {
			console->error("Can't read Squeakr file. It was written on a different endian machine.");
			exit(1);
		}
		if (cqfs.size() == 0)
			kmer_size = config.kmer_size;
		else {
			if (kmer_size != config.kmer_size) {
				console->error("Squeakr file {} has a different k-mer size. Expected: {} Available: {}",
											 squeakr_file, kmer_size, config.kmer_size);
				exit(1);
			}
		}
		if (config.cutoff == 1) {
			// console->warn("Squeakr file {} is not filtered.", squeakr_file);
		}

    cqfs.emplace_back(squeakr_file, CQF_MMAP);
		//std::string sample_id = first_part(first_part(last_part(squeakr_file, '/'),
																									//'.'), '_');
		std::string sample_id = squeakr_file;
		// console->info("Reading CQF {} Seed {}",nqf, cqfs[nqf].seed());
		// console->info("Sample id {}", sample_id);
		// cqfs.back().dump_metadata();
    inobjects.emplace_back(&cqfs[nqf], sample_id, nqf);
		if (!cqfs.front().check_similarity(&cqfs.back())) {
			console->error("Passed Squeakr files are not similar.", squeakr_file);
			exit(1);
		}
    nqf++;
	}

	ColoredDbg<SampleObject<CQF<KeyObject>*>, KeyObject> cdbg(opt.qbits,
																														inobjects[0].obj->keybits(),
																														cqfs[0].hash_mode(),
																														inobjects[0].obj->seed(),
																														prefix, nqf, MANTIS_DBG_ON_DISK);
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
								cdbg.get_cqf()->dist_elts(), cdbg.get_num_eqclasses());

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



int compare_indices_main(CompareIndicesOpt &opt)
{
	spdlog::logger* console = opt.console.get();

	console -> info("Comparing the colored dBG's.");
	
	std::string dir1(opt.cdbg1);
	if (dir1.back() != '/')
		dir1 += '/';

	std::string dir2(opt.cdbg2);
	if (dir2.back() != '/')
		dir2 += '/';


	std::string cqf1 = dir1 + mantis::CQF_FILE;
	std::vector<std::string> colClsFiles1 = mantis::fs::GetFilesExt(dir1.c_str(), mantis::EQCLASS_FILE);
	std::string sampleList1 = dir1 + mantis::SAMPLEID_FILE;
	

	console -> info("Loading the first CdBG from directory {} into memory.", dir1);
	ColoredDbg<SampleObject<CQF<KeyObject> *>, KeyObject> cdbg1(cqf1, colClsFiles1,
																sampleList1, MANTIS_DBG_IN_MEMORY);
	
	console -> info("Loaded the first CdBG; it has {} k-mers, {} color-class files, and {} color-classes.",
					cdbg1.get_cqf() -> dist_elts(), colClsFiles1.size(), cdbg1.get_num_bitvectors());


	std::string cqf2 = dir2 + mantis::CQF_FILE;
	std::vector<std::string> colClsFiles2 = mantis::fs::GetFilesExt(dir2.c_str(), mantis::EQCLASS_FILE);	  
	std::string sampleList2 = dir2 + mantis::SAMPLEID_FILE;
	
	console -> info("Loading the second CdBG from directory {} into memory.", dir2);
	ColoredDbg<SampleObject<CQF<KeyObject> *>, KeyObject> cdbg2(cqf2, colClsFiles2,
																sampleList2, MANTIS_DBG_IN_MEMORY);

	console -> info("Loaded the second CdBG; it has {} k-mers, {} color-class files, {} color-classes, and {} color-classes per bitvector buffer.",
					cdbg2.get_cqf() -> dist_elts(), colClsFiles2.size(), cdbg2.get_num_bitvectors(),
					cdbg2.get_color_class_per_buffer());


	if(cdbg1.get_cqf() -> dist_elts() != cdbg2.get_cqf() -> dist_elts())
	{
		console -> error("Mismatching number of k-mers ({}, {}).",
						cdbg1.get_cqf() -> dist_elts(), cdbg2.get_cqf() -> dist_elts());
		exit(1);
	}
	else if(cdbg1.get_num_samples() != cdbg2.get_num_samples())
	{
		console -> error("Mismatching number of samples ({}, {}).",
						cdbg1.get_num_samples(), cdbg2.get_num_samples());
		exit(1);
	}
	else
		console -> info("k-mer and sample count matches.");
																	

	auto eqclasses2 = cdbg2.get_eqclasses(), eqclasses1 = cdbg1.get_eqclasses();		
	
	// Linear scan over the CQFs
	auto it2 = cdbg2.get_cqf() -> begin(), it1 = cdbg1.get_cqf() -> begin();
	uint64_t kmerCount = 0;
	const uint64_t PROGRESS_STEP = 10000000;


	while(!it2.done() && !it1.done())
	{
		KeyObject entry2 = it2.get_cur_hash(), entry1 = it1.get_cur_hash();

		if(entry2.key != entry1.key)
		{
			console -> error("Mismatching k-mers found.");
			exit(1);
		}


		uint64_t colorId2 = entry2.count, colorId1 = entry1.count;
		std::unordered_set<std::string> sampleSet2, sampleSet1;
		const uint64_t wordLen = 64;


		uint64_t sampleCount = cdbg2.get_num_samples();

		// Fetch the set of samples in the first CdBG containing the k-mer of this iteration (entry1.key).
		uint64_t bucketIdx1 = (colorId1 - 1) / mantis::NUM_BV_BUFFER; // Assuming it's of the old mantis version.
		// uint64_t bucketIdx1 = (colorId1 - 1) / cdbg1.get_color_class_per_buffer();
		uint64_t offset1 = ((colorId1 - 1) % mantis::NUM_BV_BUFFER) * sampleCount; // Assuming it's of the old mantis version.
		// uint64_t offset1 = ((colorId1 - 1) % cdbg1.get_color_class_per_buffer()) * sampleCount;

		for(uint32_t wordCount = 0; wordCount <= sampleCount / wordLen; ++wordCount)
		{
			uint64_t readLen = std::min(wordLen, sampleCount - wordCount * wordLen);
			uint64_t word = eqclasses1[bucketIdx1].get_int(offset1, readLen);

			for(uint32_t bitIdx = 0, sampleID = wordCount * wordLen; bitIdx < readLen; ++bitIdx, ++sampleID)
				if((word >> bitIdx) & 0x01)
					sampleSet1.insert(cdbg1.get_sample(sampleID));

			offset1 += readLen;
		}
	
		// Fetch the set of samples in the second CdBG containing the k-mer of this iteration (entry2.key).
		// uint64_t bucketIdx2 = (colorId2 - 1) / mantis::NUM_BV_BUFFER;
		uint64_t bucketIdx2 = (colorId2 - 1) / cdbg2.get_color_class_per_buffer(); // Assuming it's of the new mantis version.
		// uint64_t offset2 = ((colorId2 - 1) % mantis::NUM_BV_BUFFER) * sampleCount;
		uint64_t offset2 = ((colorId2 - 1) % cdbg2.get_color_class_per_buffer()) * sampleCount; // Assuming it's of the new mantis version.

		for(uint32_t wordCount = 0; wordCount <= sampleCount / wordLen; ++wordCount)
		{
			uint64_t readLen = std::min(wordLen, sampleCount - wordCount * wordLen);
			uint64_t word = eqclasses2[bucketIdx2].get_int(offset2, readLen);

			for(uint32_t bitIdx = 0, sampleID = wordCount * wordLen; bitIdx < readLen; ++bitIdx, ++sampleID)
				if((word >> bitIdx) & 0x01)
					sampleSet2.insert(cdbg2.get_sample(sampleID));

			offset2 += readLen;
		}


		// Match the sample sets containing the k-mer of this iteration.

		if(sampleSet1.size() != sampleSet2.size())
		{
			console -> error("For one or more k-mers, the sample set sizes are found to be different between the CdBGs.");
			exit(1);
		}
		
		for(auto sample : sampleSet1)
			if(sampleSet2.find(sample) == sampleSet2.end())
			{
				console -> error("Sample sets mismatch for one or more kmers.");
				exit(1);
			}


		++it2, ++it1, kmerCount++;

		if(kmerCount % PROGRESS_STEP == 0)
			console -> info("{}M k-mers matched.", kmerCount * 10 / PROGRESS_STEP);
	}


	console -> info("CQF sizes = {}, matching k-mers found {}.", cdbg1.get_cqf() -> dist_elts(), kmerCount);
	console -> info("Matching completed. The CdBGs contain the same set of k-mers and equivalent color-classes for those k-mers.");

	return EXIT_SUCCESS;
}