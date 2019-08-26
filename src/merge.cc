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

	
	console -> info("Reading the first input colored dBG from disk.");

	std::string cqfFile1(dir1 + mantis::CQF_FILE);
	std::string sampleListFile1(dir1 + mantis::SAMPLEID_FILE);
	std::vector<std::string> colorClassFiles1 = mantis::fs::GetFilesExt(dir1.c_str(), mantis::EQCLASS_FILE);

	ColoredDbg<SampleObject<CQF<KeyObject> *>, KeyObject> cdbg1(cqfFile1, sampleListFile1, colorClassFiles1,
																MANTIS_DBG_ON_DISK);

	
	console -> info("Read colored dBG with {} k-mers and {} color-class files.",
					cdbg1.get_cqf() -> dist_elts(), cdbg1.get_eq_class_files().size());


	console -> info("Reading the second input colored dBG from disk.");

	std::string cqfFile2(dir2 + mantis::CQF_FILE);
	std::string sampleListFile2(dir2 + mantis::SAMPLEID_FILE);
	std::vector<std::string> colorClassFiles2 = mantis::fs::GetFilesExt(dir2.c_str(), mantis::EQCLASS_FILE);

	ColoredDbg<SampleObject<CQF<KeyObject> *>, KeyObject> cdbg2(cqfFile2, sampleListFile2, colorClassFiles2,
																MANTIS_DBG_ON_DISK);

	console -> info("Read colored dBG with {} k-mers and {} color-class files.",
					cdbg2.get_cqf() -> dist_elts(), cdbg2.get_eq_class_files().size());


	if(!cdbg1.get_cqf() -> check_similarity(cdbg2.get_cqf()))
	{
		console -> error("The CQF files of the colored dBGs are not similar.");
		exit(1);
	}



	ColoredDbg<SampleObject<CQF<KeyObject> *>, KeyObject> mergedCdBG(cdbg1, cdbg2, outDir, MANTIS_DBG_ON_DISK);

	if(opt.threadCount > 1)
		mergedCdBG.set_thread_count(opt.threadCount);

	if(opt.maxMemory > 1)
		mergedCdBG.set_max_memory_for_sort(opt.maxMemory);

	mergedCdBG.set_console(console);

	
	console -> info("Constructing the merged colored dBG.");

	uint64_t colorClassCount = mergedCdBG.merge_2(cdbg1, cdbg2);

	console -> info("Merged colored dBG has {} k-mers and {} color-classes.",
					mergedCdBG.get_cqf() -> dist_elts(), colorClassCount);


	mergedCdBG.get_cqf() -> dump_metadata();

	console -> info("Serializing the CQF.");

	mergedCdBG.serialize_cqf_and_abundance_dist(cdbg1, cdbg2);

	console -> info("Serialization done.");

	return EXIT_SUCCESS;
}





int validate_merge_main(ValidateMergeOpts &opt)
{
	spdlog::logger* console = opt.console.get();

	console -> info("Checking correctness of the merged CQF and color-class table.");

	
	std::string correctResPrefix(opt.correctRes);
	if (correctResPrefix.back() != '/')
		correctResPrefix += '/';

	
	std::string mergedResPrefix(opt.mergeRes);
	if (mergedResPrefix.back() != '/')
		mergedResPrefix += '/';


	std::string corrCQFfile(correctResPrefix + mantis::CQF_FILE);
	std::vector<std::string> corrEqClassFiles = mantis::fs::GetFilesExt(correctResPrefix.c_str(),
																		mantis::EQCLASS_FILE);
	std::string corrSampleList(correctResPrefix + mantis::SAMPLEID_FILE);
	

	console -> info("Loading the correct CdBG.");
	ColoredDbg<SampleObject<CQF<KeyObject>*>, KeyObject> corrCdBG(corrCQFfile, corrEqClassFiles,
																	corrSampleList, MANTIS_DBG_IN_MEMORY);
	
	console -> info("Loaded the correct CdBG; it has {} k-mers, {} color-class files, and {} color-classes.",
					corrCdBG.get_cqf() -> dist_elts(), corrEqClassFiles.size(),
					corrCdBG.get_num_bitvectors());



	std::string mergedCQFfile = mergedResPrefix + mantis::CQF_FILE;
	std::vector<std::string> mergedEqClassFiles = mantis::fs::GetFilesExt(mergedResPrefix.c_str(),
																			mantis::EQCLASS_FILE);	  
	std::string mergedSampleList = mergedResPrefix + mantis::SAMPLEID_FILE;
	
	console -> info("Loading the merged CdBG.");
	ColoredDbg<SampleObject<CQF<KeyObject>*>, KeyObject> mergedCdBG(mergedCQFfile, mergedEqClassFiles,
																	mergedSampleList, MANTIS_DBG_IN_MEMORY);

	console -> info("Loaded the merged CdBG; it has {} k-mers, {} color-class files, and {} color-classes.",
					mergedCdBG.get_cqf() -> dist_elts(), mergedEqClassFiles.size(),
					mergedCdBG.get_num_bitvectors());


	if(corrCdBG.get_cqf() -> dist_elts() != mergedCdBG.get_cqf() -> dist_elts() ||
		corrEqClassFiles.size() != mergedEqClassFiles.size() ||
		corrCdBG.get_num_bitvectors() != mergedCdBG.get_num_bitvectors())
	{
		console -> error("Mismatching meta-info.");
		exit(1);
	}
	else
		console -> info("Meta information matches.");

																	

	auto eqclass_m = mergedCdBG.get_eqclasses(), eqclass_c = corrCdBG.get_eqclasses();		


	// debug
	// for(int i = 1; i <= 3; ++i)
	// {
	// 	const uint64_t wordLen = 64;
	
	// 	uint64_t colCount_m = mergedCdBG.get_num_samples();
	// 	uint64_t bucketIdx_m = (i - 1) / mantis::NUM_BV_BUFFER;
	// 	uint64_t offset_m = ((i - 1) % mantis::NUM_BV_BUFFER) * colCount_m;

	// 	for(uint32_t wordCount = 0; wordCount <= colCount_m / wordLen; ++wordCount)
	// 	{
	// 		uint64_t readLen = std::min(wordLen, colCount_m - wordCount * wordLen);
	// 		uint64_t word = eqclass_m[bucketIdx_m].get_int(offset_m, readLen);

	// 		printf("bitvector for eq id %d is = ", (int)i);
	// 		// printf("eqId1 %d, read length = %d, word %d\n", (int)eqID_m, (int)readLen, (int)word);

	// 		// Optimize here; preferrably eliminate the loop with one statement (some sort of set_int() ?).
	// 		for(uint32_t bitIdx = 0, sampleID = wordCount * wordLen; bitIdx < readLen; ++bitIdx, ++sampleID)
	// 			if((word >> bitIdx) & 0x01)
	// 				putchar('1');
	// 			else
	// 				putchar('0');
	// 		putchar('\n');

	// 		offset_m += readLen;
	// 	}
	// }
	
	
	// Linear scan over the CQFs

	auto it_m = mergedCdBG.get_cqf() -> begin(), it_c = corrCdBG.get_cqf() -> begin();
	uint64_t kmerCount = 0;
	const uint64_t PROGRESS_STEP = 10000000;


	// For debugging purpose(s)
	uint64_t absent = 0, emptyBitVec_m = 0, emptyBitVec_c = 0;//, two = 0;


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


		// For debugging purpose(s)
		if(!eqID_m || !eqID_c)
			absent++;

		
		const uint64_t wordLen = 64;
	

		// Fetch the set of samples in the merged CdBG containing the k-mer of this iteration (entry_m.key).
		uint64_t colCount_m = mergedCdBG.get_num_samples();
		uint64_t bucketIdx_m = (eqID_m - 1) / mantis::NUM_BV_BUFFER;
		uint64_t offset_m = ((eqID_m - 1) % mantis::NUM_BV_BUFFER) * colCount_m;

		for(uint32_t wordCount = 0; wordCount <= colCount_m / wordLen; ++wordCount)
		{
			uint64_t readLen = std::min(wordLen, colCount_m - wordCount * wordLen);
			uint64_t word = eqclass_m[bucketIdx_m].get_int(offset_m, readLen);

			// printf("eqId1 %d, read length = %d, word %d\n", (int)eqID_m, (int)readLen, (int)word);

			for(uint32_t bitIdx = 0, sampleID = wordCount * wordLen; bitIdx < readLen; ++bitIdx, ++sampleID)
				if((word >> bitIdx) & 0x01)
					sampleSet_m.insert(mergedCdBG.get_sample(sampleID));

			offset_m += readLen;
		}



		// Fetch the set of samples in the correct CdBG containing the k-mer of this iteration (entry_c.key).
		uint64_t colCount_c = corrCdBG.get_num_samples();
		uint64_t bucketIdx_c = (eqID_c - 1) / mantis::NUM_BV_BUFFER;
		uint64_t offset_c = ((eqID_c - 1) % mantis::NUM_BV_BUFFER) * colCount_c;

		for(uint32_t wordCount = 0; wordCount <= colCount_c / wordLen; ++wordCount)
		{
			uint64_t readLen = std::min(wordLen, colCount_c - wordCount * wordLen);
			uint64_t word = eqclass_c[bucketIdx_c].get_int(offset_c, readLen);

			// printf("eqId2 %d, read length = %d, word %d\n", (int)eqID_c, (int)readLen, (int)word);

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

		// For debugging purpose(s)
		// if(sampleSet_m.size() == 2)
		// 	two++;

		
		// For debugging purpose(s)
		if(sampleSet_m.empty())
			emptyBitVec_m++;

		if(sampleSet_c.empty())
			emptyBitVec_c++;

		
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

	// For debugging purpose(s)
	// console -> info("Kmers found to be present in both samples = {}", two);


	console -> info("CQF sizes = {}, matching k-mers found {}.", corrCdBG.get_cqf() -> dist_elts(), kmerCount);

	// console -> info("Empty bitvectors in merged color table = {}, in correct color table = {}",
	// 					emptyBitVec_m, emptyBitVec_c);


	console -> info("Merged CdBG has correct CQF, and correct color-classes for all k-mers.");

	return EXIT_SUCCESS;
}