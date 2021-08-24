/*
 * ============================================================================
 *
 *         Author:  Jamshed Khan, jamshed@umd.edu
 *   Organization:  University of Maryland
 *
 * ============================================================================
 */

#ifndef LSMT_H
#define LSMT_H



#include "cqfMerger.h"
#include "kmer.h"


// Forward declarations.
int build_blockedCQF_main(BuildOpts& opt);
//int build_main (BuildOpts& opt);
int merge_main(MergeOpts &opt);
int build_mst_main(QueryOpts &opt);
std::vector<std::string> loadSampleFile(const std::string &sampleFileAddr);
void output_results(MSTQuery &mstQuery, std::ofstream &opfile, const std::vector<std::string> &sampleNames, QueryStats &queryStats, uint64_t numQueries);



template<typename qf_obj, typename key_obj>
class LSMT
{
    using ColoredDbg_t = ColoredDbg<qf_obj, key_obj>;
    using query_aggregate_t = std::vector<std::vector<std::pair<std::string, uint64_t>>>;
    
    public:
        // Checks if all the required data for the mantis indices at each LSMT-level
        // exists at the LSMT directory `dir`.
        static bool is_valid_LSMT(const std::string& dir, spdlog::logger* console);

        // Constructs an LSM tree configuration from the LSMT directory `dir`.
        LSMT(const std::string& dir);

        void set_console(const std::shared_ptr<spdlog::logger>& c);

        // Prints the LSM tree parameters.
        void print_config();

        // Updates the LSM tree by merging the input samples present whose paths are
        // present at the file named `sampleList`, using `threadCount` number of threads.
        void update(std::vector<std::string>& sampleList, uint threadCount);

        // Returns the k value for the k-mers that this LSM tree is built upon.
        uint64_t kmer_len();

        // Queries the k-mers from 'kmerSets' into the LSM tree;
        // where 'kmerSets' is a collection of sets of k-mers with each set corresponding
        // to one read; and writes the query results (histogram) for each read into the
        // 'output' file.
        void query(std::vector<std::unordered_set<uint64_t>> &kmerSets, std::string &output);

        void query(const QueryOpts& opt);

        void query_level(uint level, const QueryOpts& opt, query_aggregate_t& aggregate_result);

        // Queries the k-mers from 'kmerSets' into the pendling samples;
        // where 'kmerSets' is a collection of sets of k-mers with each set corresponding
        // to one read; and appends the query results (histogram) for each read into the
        // list 'resultSet'.
        void query_pending_list(std::vector<std::unordered_set<uint64_t>> &kmerSets,
                                std::vector<std::vector<std::pair<std::string, uint64_t>>> &totalResult);

        void query_pending_list(const QueryOpts& opt, query_aggregate_t& agrregate_result);


    private:
        std::string dir;
        uint scaling_factor;
        uint64_t kmer_threshold;
        uint32_t cqf_count_threshold;
        uint64_t sample_threshold;
        uint levels;
        uint sample_count = 0;
        uint qBit_init_build;
        std::vector<std::string> pending_samples;

        std::shared_ptr<spdlog::logger> shared_console;
        spdlog::logger* console;



        // Merges the Mantis index present at directory `upper_level_dir` into
        // the Mantis index (if exists) at LSM-tree[`level`]; or move the index
        // from `upper_level_dir` to LSM-tree[`level`] if currently there is no
        // Mantis index at `level`. Also, the Mantis index at `upper_level_dir`
        // is removed from disk. Uses `thread_count` number of threads in the merge.
        void merge_into_level(uint level, std::string upper_level_dir, uint thread_count);
        
        // Writes the pending list into disk.
        void dump_pending_list();

        // Builds a Mantis index for the samples collection whose paths are present
        // at the file named `input_list`, into the directory `out_dir`.
        void build_index(std::string input_list, std::string out_dir);

        // Converts the classic Mantis index present at the directory `out_dir` to
        // its MST-based version, using `thread_count` number of threads.
        void build_mst_index(std::string out_dir, uint thread_count);

        // Merges the Mantis (MST-based) indices present at the directories named
        // `cdbg1` and `cdbg2` into the directory `cdbg`, using `threadCount` number
        // of threads.
        void merge_indices(std::string cdbg1, std::string cdbg2, std::string cdbg, uint thread_count);

        // Returns the count of k-mers stored at LSM-tree[`level`].
        uint64_t kmer_count(uint level);

        // Returns the number of CQFs stored at LSM-tree[`level`].
        uint32_t cqf_count(uint level) const;

        void dump_parameters();

        bool propagation_triggered(uint64_t kmers_count, uint64_t kmers_threshold, uint32_t cqf_count, uint32_t cqf_threshold) const;

        // Remove the indices that have been merged into one.
        void remove_old_indices(std::string cdbg_dir1, std::string cdbg_dir2);

        void aggregate_query_results(query_aggregate_t& aggregate, const MSTQuery &mstQuery, const std::vector<std::string> &sampleNames, QueryStats &queryStats, uint64_t numQueries);

        void print_query_results(const query_aggregate_t& aggregate, const std::string& op_file_path) const;
};



template<typename qf_obj, typename key_obj>
bool LSMT<qf_obj, key_obj>::is_valid_LSMT(const std::string& dir, spdlog::logger* console)
{
	// Check for the LSMT parameters configuration file.
	if(!mantis::fs::FileExists((dir + mantis::PARAM_FILE).c_str()))
	{
		console->error("LSMT parameters configuration file {} does not exist in input directory {}.", mantis::PARAM_FILE, dir);
		return false;
	}

    // Check for the list of the pending samples.
	if(!mantis::fs::FileExists((dir + mantis::PENDING_SAMPLES_LIST).c_str()))
	{
		console->error("LSMT pending samples list {} does not exist in input directory {}.", mantis::PENDING_SAMPLES_LIST, dir);
		return false;
	}

    // Check validity of each LSM tree level.
	std::ifstream param_file(dir + mantis::PARAM_FILE);
    nlohmann::json params;

    param_file >> params;
    param_file.close();

	uint levels = params["levels"];

	for(uint i = 0; i < levels; ++i)
	{
		std::string level_dir = dir + mantis::LSMT_LEVEL_DIR + std::to_string(i) + "/";
        if(mantis::fs::DirExists(level_dir.c_str()) && !ColoredDbg_t::data_exists(level_dir, console))
        {
            console->error("Level {} has missing files in directory {}.", i, level_dir);
            return false;
        }
	}


	return true;
}



template<typename qf_obj, typename key_obj>
LSMT<qf_obj, key_obj>::LSMT(const std::string& dir)
{
    std::ifstream param_file(dir + mantis::PARAM_FILE);
    if(!param_file.is_open())
    {
        console->error("LSMT parameters file {} cannot be opened. Aborting.", dir + mantis::PARAM_FILE);
        std::exit(EXIT_FAILURE);
    }


    // Read in the LSM tree configurations.

    nlohmann::json param;
    param_file >> param;
    param_file.close();

    this->dir = dir,
    scaling_factor = param["scaling_factor"],
    kmer_threshold = param["kmer_threshold"],
    cqf_count_threshold = param["cqf_count_threshold"],
    sample_threshold = param["sample_threshold"],
    levels = param["levels"],
    sample_count = param["sample_count"],
    qBit_init_build = param["qBit_init_build"];


    // Load the path names of the pending samples yet to be inserted into the tree.

    std::ifstream pending_list(dir + mantis::PENDING_SAMPLES_LIST);
    if(!pending_list.is_open())
    {
        console -> error("Error reading from the pending samples list file {}. Aborting.",
                        dir + mantis::PENDING_SAMPLES_LIST);
        std::exit(EXIT_FAILURE);
    }

    std::string sample;
    while(pending_list >> sample)
        pending_samples.push_back(sample);

    pending_list.close();
}


template<typename qf_obj, typename key_obj>
void LSMT<qf_obj, key_obj>::set_console(const std::shared_ptr<spdlog::logger>& c)
{
    shared_console = c;
    console = shared_console.get();
}



template<typename qf_obj, typename key_obj>
void LSMT<qf_obj, key_obj>::print_config()
{
    console->info("LSM tree configuration:");
    console->info("LSM tree home directory: {}", dir);
    console->info("kmer threshold for level 0: {}", kmer_threshold);
    console->info("CQF count threshold for level 0: {}", cqf_count_threshold);
    console->info("Scaling factor between levels: {}", scaling_factor);
    console->info("Maximum number of pending samples at a given time: {}", sample_threshold);
    console->info("Quotient bit for the initial mantis build at an update: {}", qBit_init_build);
    console->info("Current number of levels in tree: {}", levels);
    console->info("Current number of samples in tree: {}", sample_count);
}



template<typename qf_obj, typename key_obj>
void LSMT<qf_obj, key_obj>::update(std::vector<std::string>& sample_list, uint thread_count)
{
    auto t_start = time(nullptr);


    // Read in the new samples' paths.

    for(auto p = sample_list.begin(); p != sample_list.end(); ++p)
    {
        pending_samples.push_back(*p);

        if(pending_samples.size() >= sample_threshold)
        {
            console->critical("Pushing {} samples into the LSM tree from the pending list.", pending_samples.size());

            dump_pending_list();

            build_index(dir + mantis::PENDING_SAMPLES_LIST, dir + mantis::TEMP_BUILD_IDX_DIR);
            build_mst_index(dir + mantis::TEMP_BUILD_IDX_DIR, thread_count);

            merge_into_level(0, dir + mantis::TEMP_BUILD_IDX_DIR, thread_count);

            
            uint curr_level = 0;
            uint64_t curr_kmer_threshold = kmer_threshold;
            uint32_t curr_cqf_count_threshold = cqf_count_threshold;


            // while(kmer_count(currLevel) >= currKmerThreshold)
            while(propagation_triggered(kmer_count(curr_level), kmer_threshold, cqf_count(curr_level), curr_cqf_count_threshold))
            {
                console->critical("For level {}, k-mer count is {} and CQF count is {}. k-mer and CQF count thresholds are {} and {}, respectively.",
                                    curr_level, kmer_count(curr_level), cqf_count(curr_level), kmer_threshold, curr_cqf_count_threshold);
                console->critical("Propagating the Mantis index present at this level downwards.");

                std::string levelDir = dir + mantis::LSMT_LEVEL_DIR + std::to_string(curr_level) + "/";
                merge_into_level(curr_level + 1, levelDir, thread_count);

                curr_level++;
                curr_kmer_threshold *= scaling_factor;
                curr_cqf_count_threshold *= scaling_factor;
            }

            
            if(levels < curr_level + 1)
                levels = curr_level + 1;

            sample_count += pending_samples.size();
            pending_samples.clear();

            dump_parameters();
        }
    }


    dump_pending_list();

    console -> info("Update completed. Total {} samples are kept in the pending list.", pending_samples.size());
    print_config();

    
    auto t_end = time(nullptr);
    console -> info("Time taken to update with {} samples = {} s. Total samples = {}.",
                    sample_list.size(), t_end - t_start, sample_count);
}



template<typename qf_obj, typename key_obj>
void LSMT<qf_obj, key_obj>::dump_pending_list()
{
    std::ofstream pendingList(dir + mantis::PENDING_SAMPLES_LIST);
    if(!pendingList.is_open())
    {
        console->error("Error writing to pending samples list file {}. Aborting.", dir + mantis::PENDING_SAMPLES_LIST);
        std::exit(EXIT_FAILURE);
    }

    
    for(auto p = pending_samples.begin(); p != pending_samples.end(); ++p)
        pendingList << *p << "\n";

    pendingList.close();
}



template<typename qf_obj, typename key_obj>
void LSMT<qf_obj, key_obj>::build_index(std::string input_list, std::string out_dir)
{
    BuildOpts buildOpts;

    buildOpts.flush_eqclass_dist = false;
    buildOpts.qbits = qBit_init_build;
    buildOpts.inlist = input_list;
    buildOpts.out = out_dir;
    buildOpts.console = shared_console;

    build_blockedCQF_main(buildOpts);
//    build_main(buildOpts);
}


template<typename qf_obj, typename key_obj>
void LSMT<qf_obj, key_obj>::build_mst_index(std::string out_dir, uint thread_count)
{
    QueryOpts qopt;

    qopt.prefix = out_dir;
    qopt.numThreads = thread_count;
    qopt.remove_colorClasses = true;
    qopt.console = shared_console;

    build_mst_main(qopt);
}



template<typename qf_obj, typename key_obj>
void LSMT<qf_obj, key_obj>::merge_into_level(uint level, std::string upper_level_dir, uint thread_count)
{
    std::string level_dir = dir + mantis::LSMT_LEVEL_DIR + std::to_string(level) + "/";

    // LSM tree `level` has not been reached yet, or is empty.
    if(!mantis::fs::DirExists(level_dir.c_str()))
    {
        // Make the level directory and move the mantis index files from `upper_level_dir` to there.

        ColoredDbg_t::move_index(upper_level_dir, level_dir, console);

        console->critical("Mantis index successfully moved from directory {} to {}.", upper_level_dir, level_dir);
    }
    else
    {
        // Merge the mantis index present at `upper_level_dir` into the index at directory `level_dir`.

        console->critical("Merging Mantis indices at directories {} and {} into directory {}.", level_dir, upper_level_dir, level_dir);

        merge_indices(level_dir, upper_level_dir, dir + mantis::TEMP_MERGE_IDX_DIR, thread_count);

        remove_old_indices(level_dir, upper_level_dir);
        ColoredDbg_t::move_index(dir + mantis::TEMP_MERGE_IDX_DIR, level_dir, console);

        console->critical("Mantis indices at directories {} and {} merged into directory {}.", level_dir, upper_level_dir, level_dir);
    }
}



template<typename qf_obj, typename key_obj>
void LSMT<qf_obj, key_obj>::merge_indices(std::string cdbg1, std::string cdbg2, std::string cdbg, uint threadCount)
{
    MergeOpts mergeOpts;

    mergeOpts.numThreads = threadCount;
    mergeOpts.timeLog = true;
    mergeOpts.removeIndices = true;
    mergeOpts.dir1 = cdbg1;
	mergeOpts.dir2 = cdbg2;
	mergeOpts.out = cdbg;
	mergeOpts.console = shared_console;

    merge_main(mergeOpts);
}



template<typename qf_obj, typename key_obj>
uint64_t LSMT<qf_obj, key_obj>::kmer_count(uint level)
{
    std::string level_dir = dir + mantis::LSMT_LEVEL_DIR + std::to_string(level) + "/";

    if(mantis::fs::DirExists(level_dir.c_str()))
    {
        // std::string cqfFile = levelDir + mantis::CQF_FILE;
        // CQF<key_obj> cqf(cqfFile, CQF_MMAP);

        // uint64_t kmerCount = cqf.dist_elts();
        // cqf.close();

        // return kmerCount;

        std::string index_info_file_name = level_dir + mantis::index_info_file_name;
        std::ifstream input(index_info_file_name.c_str());
        if(!input.is_open())
        {
            console->error("Error opening index info {}. Aborting.\n", index_info_file_name);
            std::exit(EXIT_FAILURE);
        }

        nlohmann::json index_info;
        input >> index_info;
        input.close();

        return index_info["num_kmers"];
    }
    
    return 0;
}



template<typename qf_obj, typename key_obj>
uint32_t LSMT<qf_obj, key_obj>::cqf_count(uint level) const
{
    std::string level_dir = dir + mantis::LSMT_LEVEL_DIR + std::to_string(level) + "/";
    std::vector<std::string> cqf_file_names = mantis::fs::GetFilesExt(level_dir.c_str(), mantis::CQF_FILE);

    return cqf_file_names.size();
}



template<typename qf_obj, typename key_obj>
uint64_t LSMT<qf_obj, key_obj>::kmer_len()
{
    for(uint level = 0; level < levels; ++level)
    {
        std::string levelDir = dir + mantis::LSMT_LEVEL_DIR + std::to_string(level) + "/";
        if(mantis::fs::DirExists(levelDir.c_str()))
        {
            std::string cqfFile = levelDir + "0_" + mantis::CQF_FILE;
            CQF<key_obj> cqf(cqfFile, CQF_FREAD);

            uint64_t kmerLen = cqf.keybits() / 2;
            cqf.close();

            return kmerLen;
        }
    }

    return 0;
}



template<typename qf_obj, typename key_obj>
void LSMT<qf_obj, key_obj>::dump_parameters()
{
    nlohmann::json paramInfo;
	{
		std::ofstream jfile(dir + mantis::PARAM_FILE);
		if (jfile.is_open())
		{
			paramInfo["dir"] = dir;
			paramInfo["scaling_factor"] = scaling_factor;
			paramInfo["kmer_threshold"] = kmer_threshold;
            paramInfo["cqf_count_threshold"] = cqf_count_threshold;
			paramInfo["sample_threshold"] = sample_threshold;
			paramInfo["levels"] = levels;
            paramInfo["sample_count"] = sample_count;
            paramInfo["qBit_init_build"] = qBit_init_build;


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
inline bool LSMT<qf_obj, key_obj>::propagation_triggered(uint64_t kmers_count, uint64_t kmers_threshold, uint32_t cqf_count, uint32_t cqf_threshold) const
{
    // return kmers_count > kmers_threshold;
    return cqf_count > cqf_threshold;
}



template<typename qf_obj, typename key_obj>
void LSMT<qf_obj, key_obj>::remove_old_indices(std::string cdbg_dir1, std::string cdbg_dir2)
{
    if(mantis::fs::DirExists(cdbg_dir1.c_str()))
    {
        std::string sys_command = std::string("rm -rf ") + cdbg_dir1;
        console->info("Removing a Mantis index from disk. System command used:");
        console->info("{}", sys_command);
        system(sys_command.c_str());
    }


    if(mantis::fs::DirExists(cdbg_dir2.c_str()))
    {
        std::string sys_command = std::string("rm -rf ") + cdbg_dir2;
        console->info("Removing a Mantis index from disk. System command used:");
        console->info("{}", sys_command);
        system(sys_command.c_str());
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
        console->error("Cannot write to file {}.", output);
        exit(1);
    }

    totalResult.resize(kmerSets.size());

    console -> info("Number of reads in query: {}.", totalResult.size());

    for(uint level = 0; level < levels; ++level)
    {
        std::string levelDir = dir + mantis::LSMT_LEVEL_DIR + std::to_string(level) + "/";
        if(mantis::fs::DirExists(levelDir.c_str()))
        {
            ColoredDbg<SampleObject<CQF<KeyObject> *>, KeyObject> cdbg(levelDir,
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
                    pending_samples.size());

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
    query(const QueryOpts& opt)
{
    auto t_start = time(nullptr);


    const std::string op_file_path(opt.output);
    std::ofstream output(op_file_path);
    if(!output.is_open())
    {
        console->error("Cannot write to file path {}.", op_file_path);
        std::exit(EXIT_FAILURE);
    }

    output.close();


    query_aggregate_t aggregate_result;


    for(uint curr_level = 0; curr_level < levels; ++curr_level)
    {
        const std::string level_dir(dir + mantis::LSMT_LEVEL_DIR + std::to_string(curr_level) + "/");
        if(mantis::fs::DirExists(level_dir.c_str()))
        {
            console->info("Querying tree level {}.", curr_level);
            query_level(curr_level, opt, aggregate_result);
        }
    }

    console->info("Query completed for full LSM-tree. Now querying the {} pending samples.", pending_samples.size());
    query_pending_list(opt, aggregate_result);
    console->info("Query completed for the pending samples.");

    console->info("Number of reads in query: {}", aggregate_result.size());


    console -> info("Serializing the query results.");
    print_query_results(aggregate_result, op_file_path);
    console -> info("Query results serialization completed.");


    auto t_end = time(nullptr);
    console -> info("Time taken to query = {} s.", t_end - t_start);
}


template<typename qf_obj, typename key_obj>
void LSMT<qf_obj, key_obj>::
    query_level(const uint level, const QueryOpts& opt, query_aggregate_t& aggregate_result)
{
    const std::string query_file(opt.query_file);
    const std::string level_dir(dir + mantis::LSMT_LEVEL_DIR + std::to_string(level) + "/");
    const std::string op_file_path(opt.output);
    
    const std::string sample_id_file_path(level_dir + mantis::SAMPLEID_FILE);
    const std::vector<std::string> sample_names = loadSampleFile(sample_id_file_path);

    QueryStats query_stats;   
    query_stats.numSamples = sample_names.size();
    console->info("Number of samples: {}", query_stats.numSamples);


    ColoredDbg<SampleObject<CQF<KeyObject> *>, KeyObject> dbg(level_dir, MANTIS_DBG_IN_MEMORY);
    dbg.set_console(console);

    console->info("Loading the first CQF.");
    const CQF<KeyObject>* const cqf = dbg.get_current_cqf();
    if(!cqf)
    {
        console->error("Failed to load CQF. Aborting.");
        std::exit(EXIT_FAILURE);
    }
    

    const uint64_t k = (cqf->keybits() / 2);
    const uint64_t query_k = (opt.k > 0 ? opt.k : k);
    console->info("Loaded the CQF. k is {}.", k);

    console->info("Loading color classes in the MST form.");
    MSTQuery mst_query(level_dir, opt.output, k, query_k, query_stats.numSamples, console, true);
    console->info("Done loading color classes. Total color class count is {}", mst_query.parentbv.size() - 1);

    console->info("Querying the colored de Bruijn graph.");
    std::ofstream output(op_file_path.c_str(), std::ofstream::app);
    
    LRUCacheMap cache_lru(100000);
    RankScores rs(1);
    uint64_t num_queries = 0;

    if(opt.process_in_bulk)
    {
        num_queries = mst_query.parseBulkKmers(query_file, k);
        console->info("# of observed sequences: {}", num_queries);
        mst_query.findSamples(dbg, cache_lru, &rs, query_stats, num_queries);

        // output_results(mst_query, output, sample_names, query_stats, num_queries);
        aggregate_query_results(aggregate_result, mst_query, sample_names, query_stats, num_queries);
    }


    output.close();
    console->info("Done querying level {}.", level);
}


template <typename qf_obj, typename key_obj>
void LSMT<qf_obj, key_obj>::
    aggregate_query_results(query_aggregate_t& aggregate, const MSTQuery &mstQuery, const std::vector<std::string> &sampleNames, QueryStats &queryStats, uint64_t numQueries)
{
    const mantis::QueryResults& result = mstQuery.allQueries;//mstQuery.getResultList(numQueries);
    if(aggregate.size() < result.size())
        aggregate.resize(result.size());

    for(auto& sample_hits: result)
    {
        // last element in the result for each query contains # of distinct kmers
        // opfile << "seq" << queryStats.cnt++ << '\t' << q[q.size()-1] << '\n';
        
        for(uint64_t i = 0; i < sample_hits.size() - 1; i++)
            if(sample_hits[i] > 0)
                aggregate[queryStats.cnt].emplace_back(sampleNames[i], sample_hits[i]);

        queryStats.cnt++;
    }
}


template <typename qf_obj, typename key_obj>
void LSMT<qf_obj, key_obj>::
    print_query_results(const query_aggregate_t& aggregate, const std::string& op_file_path) const
{
    std::ofstream output(op_file_path, std::ofstream::out);

    uint64_t seq_count = 0;
    for(auto& seq_query_result: aggregate)
    {
        output << "seq" << seq_count++ << "\n";
        for(auto& sample_hit: seq_query_result)
            output << sample_hit.first << "\t" << sample_hit.second << "\n";
    }


    output.close();
    if(output.bad())
    {
        console->error("Error writing query results.");
        std::exit(EXIT_FAILURE);
    }
}



template<typename qf_obj, typename key_obj>
void LSMT<qf_obj, key_obj>::
    query_pending_list(std::vector<std::unordered_set<uint64_t>> &kmerSets,
                        std::vector<std::vector<std::pair<std::string, uint64_t>>> &totalResult)
{
    // Go over each sample.
    for(auto sample: pending_samples)
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


template<typename qf_obj, typename key_obj>
void LSMT<qf_obj, key_obj>::
    query_pending_list(const QueryOpts& opt, query_aggregate_t& aggregate_result)
{
    const std::string query_file(opt.query_file);
    std::unordered_map<uint64_t, uint64_t> unique_kmers;
    uint64_t total_kmers = 0;
    
    std::vector<std::unordered_set<uint64_t>> kmerSets = Kmer::parse_kmers( query_file.c_str(), kmer_len(),
																			total_kmers, true,
																			unique_kmers);
    // Go over each sample.
    for(auto sample: pending_samples)
    {
        CQF<key_obj> cqf(sample, CQF_FREAD);

        // Go over each read.
        for(auto read = 0; read < kmerSets.size(); ++read)
        {
            uint64_t hitCount = 0;

            // Go over each k-mer of this read.
            for(auto kmer: kmerSets[read])
            {
                key_obj key(kmer, 0, 0);
                if(cqf.query(key, 0))
                    hitCount++;        
            }
            
            if(hitCount)
                aggregate_result[read].emplace_back(sample, hitCount);//, std::cout << "Found " << hitCount << "hits\n"; 
        }
    }
}



#endif