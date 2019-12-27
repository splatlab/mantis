/*
 * ============================================================================
 *
 *         Author:  Jamshed Khan, jamshed@cs.umd.edu
 *   Organization:  University of Maryland
 *
 * ============================================================================
 */

#ifndef _LSMT_H_
#define _LSMT_H_

#include "CdBG_Merger.h"


int build_blockedCQF_main(BuildOpts& opt);
//int build_main (BuildOpts& opt);
int merge_main(MergeOpts &opt);
int build_mst_main(QueryOpts &opt);



template<typename qf_obj, typename key_obj>
class LSMT
{
    using ColoredDbg_t = ColoredDbg<qf_obj, key_obj>;
    
    public:
        // Checks if all the required data for the mantis indices at each LSMT-level
        // exists at the LSMT directory 'dir'.
        static bool is_valid_LSMT(std::string &dir, spdlog::logger *console);

        // Constructor to load an LSM tree configuration from the LSMT directory
        // 'dir'.
        LSMT(std::string &dir);

        inline void set_console(std::shared_ptr<spdlog::logger> &c)
        { sharedConsole = c; console = sharedConsole.get(); }

        // Print the LSM tree parameters.
        void print_config();

        // Update the LSM tree by merging the input samples present, whose paths are
        // present at the list 'sampleList'.
        void update(std::vector<std::string> &sampleList, uint threadCount);

        // Returns the length of the k-mers that this LSM tree is built upon.
        uint64_t kmer_len();

        // Queries the k-mers from 'kmerSets' into the LSM tree;
        // where 'kmerSets' is a collection of sets of k-mers with each set corresponding
        // to one read; and writes the query results (histogram) for each read into the
        // 'output' file.
        void query(std::vector<std::unordered_set<uint64_t>> &kmerSets, std::string &output);

        // Queries the k-mers from 'kmerSets' into the pendling samples;
        // where 'kmerSets' is a collection of sets of k-mers with each set corresponding
        // to one read; and appends the query results (histogram) for each read into the
        // list 'resultSet'.
        void query_pending_list(std::vector<std::unordered_set<uint64_t>> &kmerSets,
                                std::vector<std::vector<std::pair<std::string, uint64_t>>> &totalResult);



    private:
        std::string dir;
        uint scalingFactor;
        uint64_t kmerThreshold;
        uint64_t sampleThreshold;
        uint levels;
        uint sampleCount = 0;
        uint qBitInitBuild;
        std::vector<std::string> pendingSamples;

        std::shared_ptr<spdlog::logger> sharedConsole;
        spdlog::logger *console;



        // Merges the mantis index present at directory 'upperLevelDir' into
        // the mantis index (if exists) at LSM-tree['level']; or move the index
        // from 'upperLevelDir' to LSM-tree['level'] if currently there is no
        // matis index at 'level'. Also, the mantis index at 'upperLevelDir'
        // is removed.
        void merge_into_level(uint level, std::string upperLevelDir, uint threadCount);
        
        void dump_pending_list();

        void build_index(std::string inputList, std::string outDir);

        void build_mst_index(std::string outDir, uint threadCount);

        void merge_indices(std::string cdbg1, std::string cdbg2, std::string cdbg, uint threadCount);

        // Returns the count of k-mers stored at LSM-tree['level'].
        uint64_t kmer_count(uint level);

        void dump_parameters();
};



template<typename qf_obj, typename key_obj>
bool LSMT<qf_obj, key_obj>::
    is_valid_LSMT(std::string &dir, spdlog::logger *console)
{
	// Check for the LSMT parameters configuration file.
	if(!mantis::fs::FileExists((dir + mantis::PARAM_FILE).c_str()))
	{
		console -> error("LSMT parameters configuration file {} does not exist in input directory {}.",
						mantis::PARAM_FILE, dir);
		return false;
	}

    // Check for the list of the pending samples.
	if(!mantis::fs::FileExists((dir + mantis::PENDING_SAMPLES_LIST).c_str()))
	{
		console -> error("LSMT pending samples list {} does not exist in input directory {}.",
						mantis::PENDING_SAMPLES_LIST, dir);
		return false;
	}

    // Check validity of each LSM tree level.
	std::ifstream paramFile(dir + mantis::PARAM_FILE);
    nlohmann::json param;

    paramFile >> param;

	uint levels = param["levels"];
	for(uint i = 0; i < levels; ++i)
	{
		std::string levelDir = dir + mantis::LSMT_LEVEL_DIR + std::to_string(i) + "/";
        if(mantis::fs::DirExists(levelDir.c_str()) && !ColoredDbg_t::data_exists(levelDir, console))
        {
            console -> error("Level {} has missing files in directory {}.", i, levelDir);
            return false;
        }
	}


	return true;
}



template<typename qf_obj, typename key_obj>
LSMT<qf_obj, key_obj>::
    LSMT(std::string &dir)
{
    std::ifstream paramFile(dir + mantis::PARAM_FILE);
    nlohmann::json param;

    paramFile >> param;
    paramFile.close();

    this -> dir = dir,
    scalingFactor = param["scalingFactor"],
    kmerThreshold = param["kmerThreshold"],
    sampleThreshold = param["sampleThreshold"],
    levels = param["levels"];
    sampleCount = param["sampleCount"];
    qBitInitBuild = param["qBitInitBuild"];


    std::ifstream pendingList(dir + mantis::PENDING_SAMPLES_LIST);
    std::string sample;

    if(!pendingList.is_open())
    {
        console -> error("Cannot read from file {}.", dir + mantis::PENDING_SAMPLES_LIST);
        exit(1);
    }

    while(pendingList >> sample)
        pendingSamples.push_back(sample);

    pendingList.close();
}



template<typename qf_obj, typename key_obj>
void LSMT<qf_obj, key_obj>::
    print_config()
{
    console -> info("LSM tree configuration:");
    console -> info("LSM tree home directory: {}", dir);
    console -> info("kmer threshold for level 0: {}", kmerThreshold);
    console -> info("Scaling factor between levels: {}", scalingFactor);
    console -> info("Maximum number of pending samples at a given time: {}", sampleThreshold);
    console -> info("Quotient bit for the initial mantis build at an update: {}", qBitInitBuild);
    console -> info("Current number of levels in tree: {}", levels);
    console -> info("Current number of samples in tree: {}", sampleCount);
}



template<typename qf_obj, typename key_obj>
void LSMT<qf_obj, key_obj>::
    update(std::vector<std::string> &sampleList, uint threadCount)
{
    auto t_start = time(nullptr);

    for(auto p = sampleList.begin(); p != sampleList.end(); ++p)
    {
        pendingSamples.push_back(*p);

        if(pendingSamples.size() >= sampleThreshold)
        {
            console -> critical("Pushing {} samples into the LSM tree from the pending list.", pendingSamples.size());

            dump_pending_list();

            build_index(dir + mantis::PENDING_SAMPLES_LIST, dir + mantis::TEMP_BUILD_IDX_DIR);
            build_mst_index(dir + mantis::TEMP_BUILD_IDX_DIR, threadCount);

            merge_into_level(0, dir + mantis::TEMP_BUILD_IDX_DIR, threadCount);

            
            uint currLevel = 0;
            uint64_t currThreshold = kmerThreshold;

            while(kmer_count(currLevel) >= currThreshold)
            {
                console -> critical("k-mer count at level {} is {}, exceeded threshold {}. Pushing this level downwards.",
                                currLevel, kmer_count(currLevel), currThreshold);

                std::string levelDir = dir + mantis::LSMT_LEVEL_DIR + std::to_string(currLevel) + "/";
                merge_into_level(currLevel + 1, levelDir, threadCount);

                currLevel++,
                currThreshold *= scalingFactor;
            }

            
            if(levels < currLevel + 1)
                levels = currLevel + 1;

            sampleCount += pendingSamples.size();
            pendingSamples.clear();

            dump_parameters();
        }
    }

    dump_pending_list();

    console -> info("Update completed. Total {} samples are kept in the pending list.", pendingSamples.size());
    print_config();

    
    auto t_end = time(nullptr);
    console -> info("Time taken to update with {} samples = {} s. Total samples = {}.",
                    sampleList.size(), t_end - t_start, sampleCount);
}



template<typename qf_obj, typename key_obj>
void LSMT<qf_obj, key_obj>::
    dump_pending_list()
{
    std::ofstream pendingList(dir + mantis::PENDING_SAMPLES_LIST);

    if(!pendingList.is_open())
    {
        console -> error("Cannot write to file {}.", dir + mantis::PENDING_SAMPLES_LIST);
        exit(1);
    }

    
    for(auto p = pendingSamples.begin(); p != pendingSamples.end(); ++p)
        pendingList << *p << "\n";

    pendingList.close();
}



template<typename qf_obj, typename key_obj>
void LSMT<qf_obj, key_obj>::
    build_index(std::string inputList, std::string outDir)
{
    BuildOpts buildOpts;

    buildOpts.flush_eqclass_dist = false;
    buildOpts.qbits = qBitInitBuild;
    buildOpts.inlist = inputList;
    buildOpts.out = outDir;
    buildOpts.console = sharedConsole;

    build_blockedCQF_main(buildOpts);
//    build_main(buildOpts);
}


template<typename qf_obj, typename key_obj>
void LSMT<qf_obj, key_obj>::
    build_mst_index(std::string outDir, uint threadCount)
{
    QueryOpts qopt;

    qopt.prefix = outDir;
    qopt.numThreads = threadCount;
    qopt.remove_colorClasses = true;
    qopt.console = sharedConsole;

    build_mst_main(qopt);
}



template<typename qf_obj, typename key_obj>
void LSMT<qf_obj, key_obj>::
    merge_into_level(uint level, std::string upperLevelDir, uint threadCount)
{
    std::string levelDir = dir + mantis::LSMT_LEVEL_DIR + std::to_string(level) + "/";

    // LSM tree 'level' does not exist, or is empty.
    if(!mantis::fs::DirExists(levelDir.c_str()))
    {
        // Make the level directory and move all the mantis index files from 'upperLevelDir' to there.

        ColoredDbg_t::move_index(upperLevelDir, levelDir, console);

        console -> critical("Mantis index successfully moved from directory {} to {}.", upperLevelDir, levelDir);
    }
    else
    {
        // Merge the mantis index present at 'upperLevelDir' into the index at directory 'levelDir'.

        merge_indices(levelDir, upperLevelDir, dir + mantis::TEMP_MERGE_IDX_DIR, threadCount);

        ColoredDbg_t::move_index(dir + mantis::TEMP_MERGE_IDX_DIR, levelDir, console);

        console -> critical("Mantis indices at directories {} and {} merged into directory {}.",
                        levelDir, upperLevelDir, levelDir);
    }
}



template<typename qf_obj, typename key_obj>
void LSMT<qf_obj, key_obj>::
    merge_indices(std::string cdbg1, std::string cdbg2, std::string cdbg, uint threadCount)
{
    MergeOpts mergeOpts;

    mergeOpts.threadCount = threadCount;
    mergeOpts.timeLog = true;
    mergeOpts.removeIndices = true;
    mergeOpts.dir1 = cdbg1;
	mergeOpts.dir2 = cdbg2;
	mergeOpts.out = cdbg;
	mergeOpts.console = sharedConsole;

    merge_main(mergeOpts);
}



template<typename qf_obj, typename key_obj>
uint64_t LSMT<qf_obj, key_obj>::
    kmer_count(uint level)
{
    std::string levelDir = dir + mantis::LSMT_LEVEL_DIR + std::to_string(level) + "/";

    if(mantis::fs::DirExists(levelDir.c_str()))
    {
        std::string cqfFile = levelDir + mantis::CQF_FILE;
        CQF<key_obj> cqf(cqfFile, CQF_MMAP);

        uint64_t kmerCount = cqf.dist_elts();
        cqf.close();

        return kmerCount;
    }
    
    return 0;
}



template<typename qf_obj, typename key_obj>
uint64_t LSMT<qf_obj, key_obj>::
    kmer_len()
{
    for(uint level = 0; level < levels; ++level)
    {
        std::string levelDir = dir + mantis::LSMT_LEVEL_DIR + std::to_string(level) + "/";
        if(mantis::fs::DirExists(levelDir.c_str()))
        {
            std::string cqfFile = levelDir + mantis::CQF_FILE;
            CQF<key_obj> cqf(cqfFile, CQF_MMAP);

            uint64_t kmerLen = cqf.keybits() / 2;
            cqf.close();

            return kmerLen;
        }
    }

    return 0;
}



template<typename qf_obj, typename key_obj>
void LSMT<qf_obj, key_obj>::
    dump_parameters()
{
    nlohmann::json paramInfo;
	{
		std::ofstream jfile(dir + mantis::PARAM_FILE);
		if (jfile.is_open())
		{
			paramInfo["dir"] = dir;
			paramInfo["scalingFactor"] = scalingFactor;
			paramInfo["kmerThreshold"] = kmerThreshold;
			paramInfo["sampleThreshold"] = sampleThreshold;
			paramInfo["levels"] = levels;
            paramInfo["sampleCount"] = sampleCount;
            paramInfo["qBitInitBuild"] = qBitInitBuild;


			jfile << paramInfo.dump(4);
		}
		else
		{
			console -> error("Could not write to configuration file {}.", dir + mantis::PARAM_FILE);
			exit(1);
		}

		console -> info("LSM-tree parameters recorded at file {}.", dir + mantis::PARAM_FILE);
		
		jfile.close();
	}
}



template<typename qf_obj, typename key_obj>
void LSMT<qf_obj, key_obj>::
    query(std::vector<std::unordered_set<uint64_t>> &kmerSets, std::string &output)
{
    auto t_start = time(nullptr);

    std::vector<std::vector<std::pair<std::string, uint64_t>>> totalResult;

    std::ofstream outputFile(output);
    if(!outputFile.is_open())
    {
        console -> error("Cannot write to file {}.", dir + mantis::PENDING_SAMPLES_LIST);
        exit(1);
    }

    totalResult.resize(kmerSets.size());

    console -> info("Number of reads in query: {}.", totalResult.size());

    for(uint level = 0; level < levels; ++level)
    {
        std::string levelDir = dir + mantis::LSMT_LEVEL_DIR + std::to_string(level) + "/";
        if(mantis::fs::DirExists(levelDir.c_str()))
        {
            std::string cqfFile = levelDir + mantis::CQF_FILE;
            std::vector<std::string> colorClassFiles = mantis::fs::GetFilesExt(levelDir.c_str(),
                                                                                mantis::EQCLASS_FILE);
            std::string sampleListFile = levelDir + mantis::SAMPLEID_FILE;
            ColoredDbg<SampleObject<CQF<KeyObject> *>, KeyObject> cdbg(cqfFile, colorClassFiles, sampleListFile,
                                                                        MANTIS_DBG_IN_MEMORY);

            uint64_t kmerLen = cdbg.get_cqf() -> keybits() / 2;

            console -> info("Loaded level {} colored dBG with {} k-mers and {} color classes.",
                            level, cdbg.get_cqf() -> dist_elts(), cdbg.get_num_bitvectors());
            
            console -> info("Querying at level {}.", level);

            // Go over each read
            for(auto k = 0; k < kmerSets.size(); ++k)
            {
                std::vector<uint64_t> result = cdbg.find_samples(kmerSets[k]);

                // Aggregate the query result for this read into this LSM tree level into
                // this read's cumulative result.
                for(auto i = 0; i < result.size(); ++i)
                    if(result[i] > 0)
                        totalResult[k].emplace_back(cdbg.get_sample(i), result[i]);
            }

            console -> info("Query done at level {}.", level);
        }
    }

    console -> info("Query completed for full LSM-tree. Now querying the {} pending samples.",
                    pendingSamples.size());

    query_pending_list(kmerSets, totalResult);

    console -> info("Query completed for the pending samples.");


    console -> info("Serializing the query results.");
    
    for(auto k = 0; k < kmerSets.size(); ++k)
    {
        outputFile << k << "\t" << kmerSets[k].size() << "\n";
        for(auto sampleCount: totalResult[k])
            outputFile << sampleCount.first << "\t" << sampleCount.second << "\n";    
    }

    console -> info("Query results serialization completed.");
    
    outputFile.flush();
    outputFile.close();

    auto t_end = time(nullptr);
    console -> info("Time taken to query with {} reads = {} s.", kmerSets.size(), t_end - t_start);
	
	// if (use_json) {
	// output_results_json(multi_kmers, cdbg, opfile, opt.process_in_bulk, uniqueKmers);
	// } else {
	// output_results(multi_kmers, cdbg, opfile, opt.process_in_bulk, uniqueKmers);
	// }
	// //std::cout << "Writing samples and abundances out." << std::endl;
	// opfile.close();
	// console->info("Writing done.");
}



template<typename qf_obj, typename key_obj>
void LSMT<qf_obj, key_obj>::
    query_pending_list(std::vector<std::unordered_set<uint64_t>> &kmerSets,
                        std::vector<std::vector<std::pair<std::string, uint64_t>>> &totalResult)
{
    // Go over each sample.
    for(auto sample: pendingSamples)
    {
        CQF<key_obj> cqf(sample, CQF_FREAD);

        // Go over each read.
        for(auto k = 0; k < kmerSets.size(); ++k)
        {
            uint64_t hitCount = 0;

            // Go over each k-mer of this read.
            for(auto kmer : kmerSets[k])
            {
                key_obj key(kmer, 0, 0);
                if(cqf.query(key, 0))
                    hitCount++;        
            }
            
            if(hitCount)
                totalResult[k].emplace_back(sample, hitCount);
        }
    }
}

#endif