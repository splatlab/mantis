/*
 * ============================================================================
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

#include "sparsepp/spp.h"
#include "tsl/sparse_map.h"

#include "MantisFS.h"
#include "ProgOpts.h"
#include "coloreddbg.h"
#include "squeakrconfig.h"
#include "json.hpp"
#include "mantis_utils.hpp"
#include "mantisconfig.hpp"

// Mantis merge: Jamshed




class MergeOpts
{
	public:
		bool flush_eqclass_dist{false};
		int qbits;
		// std::string inlist;
		std::string out;
		// int numthreads{1};
		std::shared_ptr<spdlog::logger> console{nullptr};
		std::string dir1;
		std::string dir2;

//   nlohmann::json to_json() {
//     nlohmann::json j;
//     j["dump_eqclass_dist"] = flush_eqclass_dist;
//     j["quotient_bits"] = qbits;
//     j["input_list"] = inlist;
//     j["output_dir"] = out;
//     j["num_threads"] = numthreads;
//     return j;
//   }
};





void build(BuildOpts &opt);
void test_merge_same_cdbg(BuildOpts &opt);
void test_merge_different_cdbg(BuildOpts &opt);
// void test_cqf(BuildOpts &bOpt, MergeOpts &mOpt);
void test_merge(BuildOpts &opt);
void validate_merge_result(BuildOpts &bOpt, MergeOpts &mOpt);
// For debugging purpose(s); has been required in finding a weird phenomenon: iterating over the
// CQF loaded from file produces random garbage the second time (with possible seg-faults);
// only the first iteration works.
int get_size_by_scanning(const ColoredDbg<SampleObject<CQF<KeyObject>*>, KeyObject> &cdbg);





// Mantis merge: Jamshed
/*
 * ===  FUNCTION  =============================================================
 *         Name:  main
 *  Description:  
 * ============================================================================
 */
	int
build_main ( BuildOpts& opt )
{
	// build(opt);

	// test_merge_same_cdbg(opt);

	// test_merge_different_cdbg(opt);

	test_merge(opt);

	// test_cqf(opt);

	// validate_merge_result(opt);

  return EXIT_SUCCESS;
}				/* ----------  end of function main  ---------- */




void build(BuildOpts &opt)
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
			console->warn("Squeakr file {} is not filtered.", squeakr_file);
		}

    cqfs.emplace_back(squeakr_file, CQF_MMAP);
		//std::string sample_id = first_part(first_part(last_part(squeakr_file, '/'),
																									//'.'), '_');
		std::string sample_id = squeakr_file;
		console->info("Reading CQF {} Seed {}",nqf, cqfs[nqf].seed());
		console->info("Sample id {}", sample_id);
		cqfs.back().dump_metadata();
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
}

// Mantis merge: Jamshed




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





int get_size_by_scanning(const ColoredDbg<SampleObject<CQF<KeyObject>*>, KeyObject> &cdbg)
{
	int c = 0;

	for(auto it = cdbg.get_cqf() -> begin(); !it.done(); ++it)
		c++;

	return c;
}





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





void validate_merge_result(BuildOpts &bOpt, MergeOpts &mOpt)
{
	spdlog::logger* console = mOpt.console.get();

	console -> info("Checking correctness of the merged CQF and color-class table.");

	
	std::string bOptPrefix(bOpt.out);
	if (bOptPrefix.back() != '/')
		bOptPrefix += '/';

	
	std::string mOptPrefix(mOpt.out);
	if (mOptPrefix.back() != '/')
		mOptPrefix += '/';

	std::string corrCQFfile(bOptPrefix + mantis::CQF_FILE);
	std::vector<std::string> corrEqClassFiles = mantis::fs::GetFilesExt(bOptPrefix.c_str(),
																		mantis::EQCLASS_FILE);
	std::string corrSampleList(bOptPrefix + mantis::SAMPLEID_FILE);
	

	console -> info("Loading the correct CdBG.");
	ColoredDbg<SampleObject<CQF<KeyObject>*>, KeyObject> corrCdBG(corrCQFfile, corrEqClassFiles,
																	corrSampleList, MANTIS_DBG_IN_MEMORY);
	
	console -> info("Loaded the correct CdBG; it has {} k-mers, {} color-class files, and {} color-classes (bitvectors).",
					corrCdBG.get_cqf() -> dist_elts(), corrEqClassFiles.size(),
					corrCdBG.get_num_bitvectors());



	std::string mergedCQFfile = mOptPrefix + mantis::CQF_FILE;
	std::vector<std::string> mergedEqClassFiles = mantis::fs::GetFilesExt(mOptPrefix.c_str(),
																			mantis::EQCLASS_FILE);	  
	std::string mergedSampleList = mOptPrefix + mantis::SAMPLEID_FILE;
	
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





void test_merge(BuildOpts &bOpt)
{
	MergeOpts opt;

	opt.console = bOpt.console;

	opt.qbits = bOpt.qbits;

	opt.dir1 = "in1";
	opt.dir2 = "in2";
	opt.out = "out";

	merge(opt);

	validate_merge_result(bOpt, opt);
}
// Mantis merge: Jamshed