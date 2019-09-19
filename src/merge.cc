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


bool data_exists(std::string &dir, spdlog::logger *console)
{
	if(!mantis::fs::FileExists((dir + mantis::CQF_FILE).c_str()))
	{
		console -> error("CQF file {} does not exist in input directory {}.", mantis::CQF_FILE, dir);
		return false;
	}

	if(!mantis::fs::FileExists((dir + mantis::SAMPLEID_FILE).c_str()))
	{
		console -> error("Sample-ID list file {} does not exist in input directory {}.",
						mantis::SAMPLEID_FILE, dir);
		return false;
	}

	if(mantis::fs::GetFilesExt(dir.c_str(), mantis::EQCLASS_FILE).empty())
	{
		console -> error("No equivalence-class file with extension {} exists in input directory {}.",
						mantis::EQCLASS_FILE, dir);
		return false;
	}


	return true;
}



int merge_main(MergeOpts &opt)
{
	spdlog::logger *console = opt.console.get();


	std::string dir1 = opt.dir1;
	if(dir1.back() != '/')	// Make sure it is a full directory.
		dir1 += '/';

	// Make sure if the first input directory exists.
	if(!mantis::fs::DirExists(dir1.c_str()))
	{
		console -> error("Input directory {} does not exist.", dir1);
		exit(1);
	}


	std::string dir2 = opt.dir2;
	if(dir2.back() != '/')	// Make sure it is a full directory.
		dir2 += '/';

	// Check to see if the second input directory exists.
	if(!mantis::fs::DirExists(dir2.c_str()))
	{
		console -> error("Input directory {} does not exist.", dir2);
		exit(1);
	}

	
	std::string outDir = opt.out;
	if(outDir.back() != '/')	// Make sure it is a full directory.
		outDir += '/';

	// Make the output directory if it doesn't exist.
	if(!mantis::fs::DirExists(outDir.c_str()))
		mantis::fs::MakeDir(outDir.c_str());
	
	// Check to see if the output dir exists now.
	if(!mantis::fs::DirExists(outDir.c_str()))
	{
		console -> error("Output dir {} could not be created.", outDir);
		exit(1);
	}


	// Check if all the required data exist in the input directories.
	if(!data_exists(dir1, console) || !data_exists(dir2, console))
		exit(1);

	
	console -> info("Loading metadata for the first input colored dBG from disk.");

	ColoredDbg<SampleObject<CQF<KeyObject> *>, KeyObject> cdbg1(dir1, MANTIS_DBG_ON_DISK);

	console -> info("Read colored dBG over {} samples, with {} k-mers and {} color-class files.",
					cdbg1.get_num_samples(), cdbg1.get_cqf() -> dist_elts(), cdbg1.get_eq_class_file_count());


	console -> info("Loading metadata for the second input colored dBG from disk.");

	ColoredDbg<SampleObject<CQF<KeyObject> *>, KeyObject> cdbg2(dir2, MANTIS_DBG_ON_DISK);

	console -> info("Read colored dBG over {} samples, with {} k-mers and {} color-class files.",
					cdbg2.get_num_samples(), cdbg2.get_cqf() -> dist_elts(), cdbg2.get_eq_class_file_count());


	if(!cdbg1.get_cqf() -> check_similarity(cdbg2.get_cqf()))
	{
		console -> error("The CQF files of the colored dBGs are not similar.");
		exit(1);
	}


	ColoredDbg<SampleObject<CQF<KeyObject> *>, KeyObject> mergedCdBG(cdbg1, cdbg2, outDir, MANTIS_DBG_ON_DISK);

	CdBG_Merger<SampleObject<CQF<KeyObject> *>, KeyObject> merger(cdbg1, cdbg2, mergedCdBG);
	merger.set_console(console);
	merger.set_thread_count(opt.threadCount);

	merger.merge();

	return EXIT_SUCCESS;
}



int validate_merge_main(ValidateMergeOpts &opt)
{
	spdlog::logger* console = opt.console.get();

	console -> info("Checking correctness of the merged colored dBG.");
	
	std::string correctResPrefix(opt.correctRes);
	if (correctResPrefix.back() != '/')
		correctResPrefix += '/';

	std::string mergedResPrefix(opt.mergeRes);
	if (mergedResPrefix.back() != '/')
		mergedResPrefix += '/';


	std::string corrCQFfile(correctResPrefix + mantis::CQF_FILE);
	std::vector<std::string> corrColorClassFiles = mantis::fs::GetFilesExt(correctResPrefix.c_str(),
																			mantis::EQCLASS_FILE);
	std::string corrSampleListFile(correctResPrefix + mantis::SAMPLEID_FILE);
	

	console -> info("Loading the correct CdBG into memory.");
	ColoredDbg<SampleObject<CQF<KeyObject> *>, KeyObject> corrCdBG(corrCQFfile, corrColorClassFiles,
																	corrSampleListFile, MANTIS_DBG_IN_MEMORY);
	
	console -> info("Loaded the correct CdBG; it has {} k-mers, {} color-class files, and {} color-classes.",
					corrCdBG.get_cqf() -> dist_elts(), corrColorClassFiles.size(),
					corrCdBG.get_num_bitvectors());


	std::string mergedCQFfile = mergedResPrefix + mantis::CQF_FILE;
	std::vector<std::string> mergedColorClassFiles = mantis::fs::GetFilesExt(mergedResPrefix.c_str(),
																			mantis::EQCLASS_FILE);	  
	std::string mergedSampleListFile = mergedResPrefix + mantis::SAMPLEID_FILE;
	
	console -> info("Loading the merged CdBG into memory.");
	ColoredDbg<SampleObject<CQF<KeyObject> *>, KeyObject> mergedCdBG(mergedCQFfile, mergedColorClassFiles,
																	mergedSampleListFile, MANTIS_DBG_IN_MEMORY);

	console -> info("Loaded the merged CdBG; it has {} k-mers, {} color-class files, and {} color-classes.",
					mergedCdBG.get_cqf() -> dist_elts(), mergedColorClassFiles.size(),
					mergedCdBG.get_num_bitvectors());


	if(corrCdBG.get_cqf() -> dist_elts() != mergedCdBG.get_cqf() -> dist_elts() ||
		corrColorClassFiles.size() != mergedColorClassFiles.size() ||
		corrCdBG.get_num_bitvectors() != mergedCdBG.get_num_bitvectors() ||
		corrCdBG.get_num_samples() != mergedCdBG.get_num_samples())
	{
		console -> error("Mismatching meta-info.");
		exit(1);
	}
	else
		console -> info("Meta information matches.");

																	

	auto eqclass_m = mergedCdBG.get_eqclasses(), eqclass_c = corrCdBG.get_eqclasses();		
	
	// Linear scan over the CQFs
	auto it_m = mergedCdBG.get_cqf() -> begin(), it_c = corrCdBG.get_cqf() -> begin();
	uint64_t kmerCount = 0;
	const uint64_t PROGRESS_STEP = 10000000;


	while(!it_m.done() && !it_c.done())
	{
		KeyObject entry_m = it_m.get_cur_hash(), entry_c = it_c.get_cur_hash();

		if(entry_m.key != entry_c.key)
		{
			console -> error("Mismatching k-mers found.");
			exit(1);
		}


		uint64_t eqID_m = entry_m.count, eqID_c = entry_c.count;
		std::unordered_set<std::string> sampleSet_m, sampleSet_c;
		const uint64_t wordLen = 64;


		uint64_t sampleCount = mergedCdBG.get_num_samples();
	
		// Fetch the set of samples in the merged CdBG containing the k-mer of this iteration (entry_m.key).
		uint64_t bucketIdx_m = (eqID_m - 1) / mantis::NUM_BV_BUFFER;
		uint64_t offset_m = ((eqID_m - 1) % mantis::NUM_BV_BUFFER) * sampleCount;

		for(uint32_t wordCount = 0; wordCount <= sampleCount / wordLen; ++wordCount)
		{
			uint64_t readLen = std::min(wordLen, sampleCount - wordCount * wordLen);
			uint64_t word = eqclass_m[bucketIdx_m].get_int(offset_m, readLen);

			for(uint32_t bitIdx = 0, sampleID = wordCount * wordLen; bitIdx < readLen; ++bitIdx, ++sampleID)
				if((word >> bitIdx) & 0x01)
					sampleSet_m.insert(mergedCdBG.get_sample(sampleID));

			offset_m += readLen;
		}


		// Fetch the set of samples in the correct CdBG containing the k-mer of this iteration (entry_c.key).
		uint64_t bucketIdx_c = (eqID_c - 1) / mantis::NUM_BV_BUFFER;
		uint64_t offset_c = ((eqID_c - 1) % mantis::NUM_BV_BUFFER) * sampleCount;

		for(uint32_t wordCount = 0; wordCount <= sampleCount / wordLen; ++wordCount)
		{
			uint64_t readLen = std::min(wordLen, sampleCount - wordCount * wordLen);
			uint64_t word = eqclass_c[bucketIdx_c].get_int(offset_c, readLen);

			for(uint32_t bitIdx = 0, sampleID = wordCount * wordLen; bitIdx < readLen; ++bitIdx, ++sampleID)
				if((word >> bitIdx) & 0x01)
					sampleSet_c.insert(corrCdBG.get_sample(sampleID));

			offset_c += readLen;
		}

		if(sampleSet_c.size() != sampleSet_m.size())
		{
			console -> error("For one or more k-mers, the sample set sizes are found to be different between the CdBGs.");
			exit(1);
		}

		
		for(auto sample : sampleSet_c)
			if(sampleSet_m.find(sample) == sampleSet_m.end())
			{
				console -> error("Sample sets mismatch for one or more kmers.");
				exit(1);
			}


		++it_m, ++it_c, kmerCount++;

		if(kmerCount % PROGRESS_STEP == 0)
			console -> info("{}M k-mers matched.", kmerCount * 10 / PROGRESS_STEP);
	}


	console -> info("CQF sizes = {}, matching k-mers found {}.", corrCdBG.get_cqf() -> dist_elts(), kmerCount);
	console -> info("Merged CdBG has correct CQF, and correct color-classes for all k-mers.");

	return EXIT_SUCCESS;
}