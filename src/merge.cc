/*
 * ============================================================================
 *
 *         Author:  Jamshed Khan (mdkhan@cs.stonybrook.edu)
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

#include "sparsepp/spp.h"
#include "tsl/sparse_map.h"

#include "MantisFS.h"
#include "ProgOpts.h"
#include "coloreddbg.h"
#include "squeakrconfig.h"
#include "json.hpp"
// #include "mantis_utils.hpp"
#include "mantisconfig.hpp"



/*
void test_merge_same_cdbg(BuildOpts &opt);
void test_merge_different_cdbg(BuildOpts &opt);
// void test_cqf(BuildOpts &bOpt, MergeOpts &mOpt);
void test_merge(BuildOpts &opt);
void validate_merge_result(BuildOpts &bOpt, MergeOpts &mOpt);
// For debugging purpose(s); has been required in finding a weird phenomenon: iterating over the
// CQF loaded from file produces random garbage the second time (with possible seg-faults);
// only the first iteration works.
int get_size_by_scanning(const ColoredDbg<SampleObject<CQF<KeyObject>*>, KeyObject> &cdbg);
*/



/*
void test_merge_same_cdbg(BuildOpts& opt)
{
	puts("\n\n========================================\nMSG: In test function.\n\n");


	std::string prefix(opt.out);
	std::string dbgFile(prefix + mantis::CQF_FILE);
	std::string sampleListFile(prefix + mantis::SAMPLEID_FILE);
	std::vector<std::string> eqclassFiles = mantis::fs::GetFilesExt(prefix.c_str(), mantis::EQCLASS_FILE);

	ColoredDbg<SampleObject<CQF<KeyObject>*>, KeyObject> inputCDBG(dbgFile, sampleListFile, eqclassFiles,
																	MANTIS_DBG_ON_DISK);

	// This inputCDBG malfunctions at the second and onward iterations of its CQF.
	// ColoredDbg<SampleObject<CQF<KeyObject>*>, KeyObject> inputCDBG(dbgFile, sampleListFile, eqclassFiles,
	// 																MANTIS_DBG_IN_MEMORY);


	// ColoredDbg<SampleObject<CQF<KeyObject>*>, KeyObject> newDBG(tempDBG, tempDBG, opt.qbits, prefix,
	// 																MANTIS_DBG_IN_MEMORY);

	char OUTPUT_CQF_FILE[] = "merged_dbg_cqf.ser";
	ColoredDbg<SampleObject<CQF<KeyObject>*>, KeyObject> mergedCdBG(inputCDBG, inputCDBG, opt.qbits, prefix,
																	MANTIS_DBG_ON_DISK);

	
	// For debugging purpose(s).
	printf("\nSize of loaded CQF: through scanning %d, through CQF size call %d\n\n",
		get_size_by_scanning(inputCDBG), (int)inputCDBG.get_cqf() -> dist_elts());
	
	//mergedCdBG.construct(sampleDBG, sampleDBG);
	mergedCdBG.construct(inputCDBG, inputCDBG);

	// For debugging purpose(s); size through iteration gets corrupted if MANTIS_DBG_IN_MEMORY is used in
	// the constructor of an input CDBG.
	printf("\nSize of loaded CQF: through scanning %d, through CQF size call %d\n\n",
		get_size_by_scanning(inputCDBG), (int)inputCDBG.get_cqf() -> dist_elts());

	
	puts("\n\nMSG: Exiting test function.\n========================================\n\n");
}
*/


/*
void test_merge_different_cdbg(BuildOpts& opt)
{
	puts("\n\n========================================\nMSG: In test function.\n\n");


	std::string prefix(opt.out);

	std::string dbgFile0(prefix + "dbg_cqf0.ser");
	std::string sampleListFile0(prefix + "sampleid0.lst");
	std::vector<std::string> eqclassFiles0 = mantis::fs::GetFilesExt(prefix.c_str(), "eqclass_rrr0.cls");

	std::string dbgFile1(prefix + "dbg_cqf1.ser");
	std::string sampleListFile1(prefix + "sampleid1.lst");
	std::vector<std::string> eqclassFiles1 = mantis::fs::GetFilesExt(prefix.c_str(), "eqclass_rrr1.cls");

	ColoredDbg<SampleObject<CQF<KeyObject>*>, KeyObject> inputCDBG0(dbgFile0, sampleListFile0, eqclassFiles0,
																	MANTIS_DBG_ON_DISK);

	printf("\nLoaded CQF0 size = %llu\n", (unsigned long long)inputCDBG0.get_cqf() -> dist_elts());

	ColoredDbg<SampleObject<CQF<KeyObject>*>, KeyObject> inputCDBG1(dbgFile1, sampleListFile1, eqclassFiles1,
																	MANTIS_DBG_ON_DISK);

	printf("\nLoaded CQF1 size = %llu\n", (unsigned long long)inputCDBG1.get_cqf() -> dist_elts());

	// This inputCDBG malfunctions at the second and onward iterations of its CQF.
	// ColoredDbg<SampleObject<CQF<KeyObject>*>, KeyObject> inputCDBG(dbgFile, sampleListFile, eqclassFiles,
	// 																MANTIS_DBG_IN_MEMORY);


	// ColoredDbg<SampleObject<CQF<KeyObject>*>, KeyObject> newDBG(tempDBG, tempDBG, opt.qbits, prefix,
	// 																MANTIS_DBG_IN_MEMORY);

	char OUTPUT_CQF_FILE[] = "merged_dbg_cqf.ser";
	ColoredDbg<SampleObject<CQF<KeyObject>*>, KeyObject> mergedCdBG(inputCDBG0, inputCDBG1, opt.qbits, prefix,
																	MANTIS_DBG_ON_DISK);

	
	// For debugging purpose(s).
	// printf("\nSize of loaded CQF: through scanning %d, through CQF size call %d\n\n",
	// 	get_size_by_scanning(inputCDBG), (int)inputCDBG.get_cqf() -> dist_elts());
	
	//mergedCdBG.construct(sampleDBG, sampleDBG);
	mergedCdBG.construct(inputCDBG0, inputCDBG1);

	mergedCdBG.serialize(inputCDBG0, inputCDBG1);

	// For debugging purpose(s); size through iteration gets corrupted if MANTIS_DBG_IN_MEMORY is used in
	// the constructor of an input CDBG.
	// printf("\nSize of loaded CQF: through scanning %d, through CQF size call %d\n\n",
	// 	get_size_by_scanning(inputCDBG), (int)inputCDBG.get_cqf() -> dist_elts());

	// if (dbg_alloc_flag == MANTIS_DBG_IN_MEMORY)
	// 	dbg.serialize(prefix + mantis::CQF_FILE);
	// else
	// 	dbg.close();

	
	puts("\n\nMSG: Exiting test function.\n========================================\n\n");
}
*/




/*
void test_cqf(BuildOpts &opt)
{
	std::string prefix = opt.out;
	std::string fileName = prefix + "merged_dbg_cqf.ser";
	CQF<KeyObject> mergedCQF(fileName, CQF_MMAP);

	printf("\nSize of loaded CDBG (merged) = %llu\n", (unsigned long long)mergedCQF.dist_elts());


	int c = 0;
	for(auto it = mergedCQF.begin(); !it.done(); ++it)
		c++;

	printf("\nSize of loaded CDBG (merged) through scanning = %d\n", c);


	fileName = prefix + mantis::CQF_FILE;
	CQF<KeyObject> correctCQF(fileName, CQF_MMAP);


	auto it_m = mergedCQF.begin(), it_c = correctCQF.begin();
	while(!it_m.done() && !it_c.done())
	{
		KeyObject mergedEntry = it_m.get_cur_hash(), correctEntry = it_c.get_cur_hash();

		if(mergedEntry.key != correctEntry.key)
		{
			puts("\nMSG: Incorrect kmer entry found.\n\n");
			break;
		}
		else
			++it_m, ++it_c;
	}

	if(it_m.done() && it_c.done())
		puts("\nMSG: All kmers match in the CQFs.\n");
	else
		puts("\nMSG: Mismatching kmer entries exist in the CQFs.\n");


	// c = 0;
	// for(auto it = cqf.begin(); !it.done(); ++it)
	// 	c++;

	// // Second time, to check if that weird phenomenon happens again.
	// printf("\nSize of loaded CDBG (merged) through scanning = %d\n", c);
}
*/





/*
int get_size_by_scanning(const ColoredDbg<SampleObject<CQF<KeyObject>*>, KeyObject> &cdbg)
{
	int c = 0;

	for(auto it = cdbg.get_cqf() -> begin(); !it.done(); ++it)
		c++;

	return c;
}
*/





bool data_exists(std::string &dir, spdlog::logger *console)
{
	if(!mantis::fs::FileExists((dir + mantis::CQF_FILE).c_str()))
	{
		console -> error("CQF file {} does not exist in input directory {}.", mantis::CQF_FILE, dir);
		return false;
	}

	if(!mantis::fs::FileExists((dir + mantis::SAMPLEID_FILE).c_str()))
	{
		console -> error("Sample ID list file {} does not exist in input directory {}.",
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





/*
void merge(MergeOpts &opt)
{
	spdlog::logger* console = opt.console.get();


	std::string dir1 = opt.dir1;
	if(dir1.back() != '/')	// Make sure it is a full directory
		dir1 += '/';

	// Make sure if the first input directory exists.
	if(!mantis::fs::DirExists(dir1.c_str()))
	{
		console -> error("Input directory {} does not exist.", dir1);
		exit(1);
	}


	std::string dir2 = opt.dir2;
	if(dir2.back() != '/')	// Make sure it is a full directory
		dir2 += '/';

	// Check to see if the second input directory exists.
	if(!mantis::fs::DirExists(dir2.c_str()))
	{
		console -> error("Input directory {} does not exist.", dir2);
		exit(1);
	}

	
	std::string outDir = opt.out;
	if(outDir.back() != '/')	// Make sure it is a full directory
		outDir += '/';

	// Make the output directory if it doesn't exist.
	if(!mantis::fs::DirExists(outDir.c_str()))
		mantis::fs::MakeDir(outDir.c_str());
	
	// Check to see if the output dir exists now.
	if(!mantis::fs::DirExists(outDir.c_str()))
	{
		console->error("Output dir {} could not be successfully created.", outDir);
		exit(1);
	}


	// Check if all the required data exist in the input directories.
	if(!data_exists(dir1, console) || !data_exists(dir2, console))
		exit(1);


	// If we made it this far, record relevant meta information in the output directory
//   nlohmann::json minfo;
//   {
//     std::ofstream jfile(prefix + "/" + mantis::meta_file_name);
//     if (jfile.is_open()) {
//       minfo = opt.to_json();
//       minfo["start_time"] = mantis::get_current_time_as_string();
//       minfo["mantis_version"] = mantis::version;
//       minfo["index_version"] = mantis::index_version;
//       jfile << minfo.dump(4);
//     } else {
//       console->error("Could not write to output directory {}", prefix);
//       exit(1);
//     }
//     jfile.close();
//   }

	
	console -> info("Reading the first input Colored dBG from disk.");

	std::string dbgFile1(dir1 + mantis::CQF_FILE);
	std::string sampleListFile1(dir1 + mantis::SAMPLEID_FILE);
	std::vector<std::string> eqclassFiles1 = mantis::fs::GetFilesExt(dir1.c_str(), mantis::EQCLASS_FILE);

	ColoredDbg<SampleObject<CQF<KeyObject>*>, KeyObject> cdbg1(dbgFile1, sampleListFile1, eqclassFiles1,
																	MANTIS_DBG_ON_DISK);

	// printf("\nCQF size = %d\n\n", (int)cdbg1.get_cqf() -> dist_elts());
	
	console -> info("Read colored dBG with {} k-mers and {} color-class files.",
					cdbg1.get_cqf() -> dist_elts(), cdbg1.get_eq_class_files().size());


	console -> info("Reading the second input Colored dBG from disk.");
	std::string dbgFile2(dir2 + mantis::CQF_FILE);
	std::string sampleListFile2(dir2 + mantis::SAMPLEID_FILE);
	std::vector<std::string> eqclassFiles2 = mantis::fs::GetFilesExt(dir2.c_str(), mantis::EQCLASS_FILE);


	ColoredDbg<SampleObject<CQF<KeyObject>*>, KeyObject> cdbg2(dbgFile2, sampleListFile2, eqclassFiles2,
																	MANTIS_DBG_ON_DISK);

	console -> info("Read colored dBG with {} k-mers and {} color-class files.",
					cdbg2.get_cqf() -> dist_elts(), cdbg2.get_eq_class_files().size());



	ColoredDbg<SampleObject<CQF<KeyObject>*>, KeyObject> mergedCdBG(cdbg1, cdbg2, opt.qbits, outDir,
																	MANTIS_DBG_ON_DISK);

	
	if(opt.flush_eqclass_dist)
		mergedCdBG.set_flush_eqclass_dist();

	
	console -> info("Constructing the merged Colored dBG.");

	mergedCdBG.construct(cdbg1, cdbg2);

	console -> info("Merged colored dBG has {} k-mers and {} equivalence classes",
								mergedCdBG.get_cqf() -> dist_elts(), mergedCdBG.get_eqclass_count());

	mergedCdBG.get_cqf() -> dump_metadata();


	console -> info("Serializing the CQF and the eq classes in directory {}.", outDir);

	mergedCdBG.serialize(cdbg1, cdbg2);

	console -> info("Serialization done.");


	//   {
//     std::ofstream jfile(prefix + "/" + mantis::meta_file_name);
//     if (jfile.is_open()) {
//       minfo["end_time"] = mantis::get_current_time_as_string();
//       jfile << minfo.dump(4);
//     } else {
//       console->error("Could not write to output directory {}", prefix);
//     }
//     jfile.close();
//   }
// }
}
*/




/*
void validate_merge_result(BuildOpts &bOpt, MergeOpts &mOpt)
{
	spdlog::logger* console = mOpt.console.get();

	console -> info("Checking correctness of the merged CQF and color-class table.");

	
	std::string correctResPrefix(bOpt.out);
	if (correctResPrefix.back() != '/')
		correctResPrefix += '/';

	
	std::string mergedResPrefix(mOpt.out);
	if (mergedResPrefix.back() != '/')
		mergedResPrefix += '/';

	std::string corrCQFfile(correctResPrefix + mantis::CQF_FILE);
	std::vector<std::string> corrEqClassFiles = mantis::fs::GetFilesExt(correctResPrefix.c_str(),
																		mantis::EQCLASS_FILE);
	std::string corrSampleList(correctResPrefix + mantis::SAMPLEID_FILE);
	

	console -> info("Loading the correct CdBG.");
	ColoredDbg<SampleObject<CQF<KeyObject>*>, KeyObject> corrCdBG(corrCQFfile, corrEqClassFiles,
																	corrSampleList, MANTIS_DBG_IN_MEMORY);
	
	console -> info("Loaded the correct CdBG; it has {} k-mers, {} color-class files, and {} color-classes (bitvectors).",
					corrCdBG.get_cqf() -> dist_elts(), corrEqClassFiles.size(),
					corrCdBG.get_num_bitvectors());



	std::string mergedCQFfile = mergedResPrefix + mantis::CQF_FILE;
	std::vector<std::string> mergedEqClassFiles = mantis::fs::GetFilesExt(mergedResPrefix.c_str(),
																			mantis::EQCLASS_FILE);	  
	std::string mergedSampleList = mergedResPrefix + mantis::SAMPLEID_FILE;
	
	console -> info("Loading the merged CdBG.");
	ColoredDbg<SampleObject<CQF<KeyObject>*>, KeyObject> mergedCdBG(mergedCQFfile, mergedEqClassFiles,
																	mergedSampleList, MANTIS_DBG_IN_MEMORY);

	console -> info("Loaded the merged CdBG; it has {} k-mers, {} color-class files, and {} color-classes(bitvectors).",
					mergedCdBG.get_cqf() -> dist_elts(), mergedEqClassFiles.size(),
					mergedCdBG.get_num_bitvectors());

																	

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


	// For debugging purpose(s)
	uint64_t c = 0, absent = 0, emptyBitVec_m = 0, emptyBitVec_c = 0;//, two = 0;


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


		++it_m, ++it_c, c++;
	}

	// For debugging purpose(s)
	// console -> info("Kmers found to be present in both samples = {}", two);


	console -> info("CQF iteration count = {}, 0-count kmers found = {}.", c, absent);

	console -> info("Empty bitvectors in merged color table = {}, in correct color table = {}",
						emptyBitVec_m, emptyBitVec_c);


	console -> info("Merged CdBG has correct CQF, and correct color classes for all k-mers.");
}
*/




/*
void test_merge(BuildOpts &bOpt)
{
	MergeOpts opt;

	opt.console = bOpt.console;
	opt.qbits = bOpt.qbits;
	opt.flush_eqclass_dist = bOpt.flush_eqclass_dist;

	opt.dir1 = "in1";
	opt.dir2 = "in2";
	opt.out = "out";

	merge(opt);

	validate_merge_result(bOpt, opt);
}
*/





int merge_main(MergeOpts &opt)
{
	spdlog::logger* console = opt.console.get();


	std::string dir1 = opt.dir1;
	if(dir1.back() != '/')	// Make sure it is a full directory
		dir1 += '/';

	// Make sure if the first input directory exists.
	if(!mantis::fs::DirExists(dir1.c_str()))
	{
		console -> error("Input directory {} does not exist.", dir1);
		exit(1);
	}


	std::string dir2 = opt.dir2;
	if(dir2.back() != '/')	// Make sure it is a full directory
		dir2 += '/';

	// Check to see if the second input directory exists.
	if(!mantis::fs::DirExists(dir2.c_str()))
	{
		console -> error("Input directory {} does not exist.", dir2);
		exit(1);
	}

	
	std::string outDir = opt.out;
	if(outDir.back() != '/')	// Make sure it is a full directory
		outDir += '/';

	// Make the output directory if it doesn't exist.
	if(!mantis::fs::DirExists(outDir.c_str()))
		mantis::fs::MakeDir(outDir.c_str());
	
	// Check to see if the output dir exists now.
	if(!mantis::fs::DirExists(outDir.c_str()))
	{
		console->error("Output dir {} could not be successfully created.", outDir);
		exit(1);
	}


	// Check if all the required data exist in the input directories.
	if(!data_exists(dir1, console) || !data_exists(dir2, console))
		exit(1);


	// If we made it this far, record relevant meta information in the output directory
//   nlohmann::json minfo;
//   {
//     std::ofstream jfile(prefix + "/" + mantis::meta_file_name);
//     if (jfile.is_open()) {
//       minfo = opt.to_json();
//       minfo["start_time"] = mantis::get_current_time_as_string();
//       minfo["mantis_version"] = mantis::version;
//       minfo["index_version"] = mantis::index_version;
//       jfile << minfo.dump(4);
//     } else {
//       console->error("Could not write to output directory {}", prefix);
//       exit(1);
//     }
//     jfile.close();
//   }

	
	console -> info("Reading the first input Colored dBG from disk.");

	std::string dbgFile1(dir1 + mantis::CQF_FILE);
	std::string sampleListFile1(dir1 + mantis::SAMPLEID_FILE);
	std::vector<std::string> eqclassFiles1 = mantis::fs::GetFilesExt(dir1.c_str(), mantis::EQCLASS_FILE);

	ColoredDbg<SampleObject<CQF<KeyObject>*>, KeyObject> cdbg1(dbgFile1, sampleListFile1, eqclassFiles1,
																	MANTIS_DBG_ON_DISK);

	
	console -> info("Read colored dBG with {} k-mers and {} color-class files.",
					cdbg1.get_cqf() -> dist_elts(), cdbg1.get_eq_class_files().size());


	console -> info("Reading the second input Colored dBG from disk.");
	std::string dbgFile2(dir2 + mantis::CQF_FILE);
	std::string sampleListFile2(dir2 + mantis::SAMPLEID_FILE);
	std::vector<std::string> eqclassFiles2 = mantis::fs::GetFilesExt(dir2.c_str(), mantis::EQCLASS_FILE);


	ColoredDbg<SampleObject<CQF<KeyObject>*>, KeyObject> cdbg2(dbgFile2, sampleListFile2, eqclassFiles2,
																	MANTIS_DBG_ON_DISK);

	console -> info("Read colored dBG with {} k-mers and {} color-class files.",
					cdbg2.get_cqf() -> dist_elts(), cdbg2.get_eq_class_files().size());


	if(!cdbg1.get_cqf() -> check_similarity(cdbg2.get_cqf()))
	{
		console -> error("The CQF files of the colored dBGs are not similar.");
		exit(1);
	}



	ColoredDbg<SampleObject<CQF<KeyObject>*>, KeyObject> mergedCdBG(cdbg1, cdbg2, outDir,
																	MANTIS_DBG_ON_DISK);

	
	if(opt.flush_eqclass_dist)
		mergedCdBG.set_flush_eqclass_dist();

	mergedCdBG.set_console(console);

	
	console -> info("Constructing the merged Colored dBG.");

	mergedCdBG.construct(cdbg1, cdbg2);

	console -> info("Merged colored dBG has {} k-mers and {} equivalence classes",
								mergedCdBG.get_cqf() -> dist_elts(), mergedCdBG.get_eqclass_count());

	mergedCdBG.get_cqf() -> dump_metadata();


	console -> info("Serializing the CQF and the eq classes in directory {}.", outDir);

	mergedCdBG.serialize(cdbg1, cdbg2);

	console -> info("Serialization done.");


	//   {
//     std::ofstream jfile(prefix + "/" + mantis::meta_file_name);
//     if (jfile.is_open()) {
//       minfo["end_time"] = mantis::get_current_time_as_string();
//       jfile << minfo.dump(4);
//     } else {
//       console->error("Could not write to output directory {}", prefix);
//     }
//     jfile.close();
//   }
// }

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





int test_merge_main(MergeOpts &opt)
{
	spdlog::logger* console = opt.console.get();


	std::string dir1 = opt.dir1;
	if(dir1.back() != '/')	// Make sure it is a full directory
		dir1 += '/';

	// Make sure if the first input directory exists.
	if(!mantis::fs::DirExists(dir1.c_str()))
	{
		console -> error("Input directory {} does not exist.", dir1);
		exit(1);
	}


	std::string dir2 = opt.dir2;
	if(dir2.back() != '/')	// Make sure it is a full directory
		dir2 += '/';

	// Check to see if the second input directory exists.
	if(!mantis::fs::DirExists(dir2.c_str()))
	{
		console -> error("Input directory {} does not exist.", dir2);
		exit(1);
	}

	
	std::string outDir = opt.out;
	if(outDir.back() != '/')	// Make sure it is a full directory
		outDir += '/';

	// Make the output directory if it doesn't exist.
	if(!mantis::fs::DirExists(outDir.c_str()))
		mantis::fs::MakeDir(outDir.c_str());
	
	// Check to see if the output dir exists now.
	if(!mantis::fs::DirExists(outDir.c_str()))
	{
		console->error("Output dir {} could not be successfully created.", outDir);
		exit(1);
	}


	// Check if all the required data exist in the input directories.
	if(!data_exists(dir1, console) || !data_exists(dir2, console))
		exit(1);


	// If we made it this far, record relevant meta information in the output directory
//   nlohmann::json minfo;
//   {
//     std::ofstream jfile(prefix + "/" + mantis::meta_file_name);
//     if (jfile.is_open()) {
//       minfo = opt.to_json();
//       minfo["start_time"] = mantis::get_current_time_as_string();
//       minfo["mantis_version"] = mantis::version;
//       minfo["index_version"] = mantis::index_version;
//       jfile << minfo.dump(4);
//     } else {
//       console->error("Could not write to output directory {}", prefix);
//       exit(1);
//     }
//     jfile.close();
//   }

	
	console -> info("Reading the first input Colored dBG from disk.");

	std::string dbgFile1(dir1 + mantis::CQF_FILE);
	std::string sampleListFile1(dir1 + mantis::SAMPLEID_FILE);
	std::vector<std::string> eqclassFiles1 = mantis::fs::GetFilesExt(dir1.c_str(), mantis::EQCLASS_FILE);

	ColoredDbg<SampleObject<CQF<KeyObject>*>, KeyObject> cdbg1(dbgFile1, sampleListFile1, eqclassFiles1,
																	MANTIS_DBG_ON_DISK);

	
	console -> info("Read colored dBG with {} k-mers and {} color-class files.",
					cdbg1.get_cqf() -> dist_elts(), cdbg1.get_eq_class_files().size());


	console -> info("Reading the second input Colored dBG from disk.");
	std::string dbgFile2(dir2 + mantis::CQF_FILE);
	std::string sampleListFile2(dir2 + mantis::SAMPLEID_FILE);
	std::vector<std::string> eqclassFiles2 = mantis::fs::GetFilesExt(dir2.c_str(), mantis::EQCLASS_FILE);


	ColoredDbg<SampleObject<CQF<KeyObject>*>, KeyObject> cdbg2(dbgFile2, sampleListFile2, eqclassFiles2,
																	MANTIS_DBG_ON_DISK);

	console -> info("Read colored dBG with {} k-mers and {} color-class files.",
					cdbg2.get_cqf() -> dist_elts(), cdbg2.get_eq_class_files().size());


	if(!cdbg1.get_cqf() -> check_similarity(cdbg2.get_cqf()))
	{
		console -> error("The CQF files of the colored dBGs are not similar.");
		exit(1);
	}



	ColoredDbg<SampleObject<CQF<KeyObject>*>, KeyObject> mergedCdBG(cdbg1, cdbg2, outDir,
																	MANTIS_DBG_ON_DISK);

	
	if(opt.flush_eqclass_dist)
		mergedCdBG.set_flush_eqclass_dist();

	if(opt.threadCount > 1)
		mergedCdBG.set_thread_count(opt.threadCount);

	if(opt.maxMemory > 1)
		mergedCdBG.set_max_memory_for_sort(opt.maxMemory);

	mergedCdBG.set_console(console);

	
	console -> info("Constructing the merged Colored dBG.");

	mergedCdBG.merge(cdbg1, cdbg2);

	console -> info("Merged colored dBG has {} k-mers and {} equivalence classes",
								mergedCdBG.get_cqf() -> dist_elts(), mergedCdBG.get_color_class_count());

	mergedCdBG.get_cqf() -> dump_metadata();


	console -> info("Serializing the CQF in directory {}.", outDir);

	// mergedCdBG.serialize(cdbg1, cdbg2);
	mergedCdBG.serialize_cqf_and_abundance_dist(cdbg1, cdbg2);

	console -> info("Serialization done.");


	//   {
//     std::ofstream jfile(prefix + "/" + mantis::meta_file_name);
//     if (jfile.is_open()) {
//       minfo["end_time"] = mantis::get_current_time_as_string();
//       jfile << minfo.dump(4);
//     } else {
//       console->error("Could not write to output directory {}", prefix);
//     }
//     jfile.close();
//   }
// }

	return EXIT_SUCCESS;
}
// Mantis merge: Jamshed