/*
 * ============================================================================
 *
 *       Filename:  dbgccmst.cc
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  2017-10-27 12:56:50 AM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Rob Johnson <robj@vmware.com>
 *   Organization:  VMware, Inc.
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
#include <unordered_map>
#include <set>
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
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/kruskal_min_spanning_tree.hpp>
#include <boost/range/adaptor/map.hpp>

#include "MantisFS.h"
#include "ProgOpts.h"
#include "spdlog/spdlog.h"
#include "kmer.h"
#include "coloreddbg.h"
#include "common_types.h"
#include "CLI/CLI.hpp"
#include "CLI/Timer.hpp"
#include "hashutil.h"

static inline int popcnt(uint64_t val)
{
	asm("popcnt %[val], %[val]"
			: [val] "+r" (val)
			:
			: "cc");
	return val;
}

typedef ColoredDbg<SampleObject<CQF<KeyObject>*>, KeyObject> colored_dbg;

using namespace boost;
typedef adjacency_list < vecS, vecS, undirectedS, no_property, property < edge_weight_t, uint64_t > > ccGraph;

// Check whether kmer2 exists and, if it does, add the corresponding
// edge to the color class graph.
//
// Assumes: kmer1 exists with color id ccid1.
static void add_edge(colored_dbg &cdbg,
										 ccGraph &ccg,
										 uint64_t ccid1,
										 uint64_t kmer2) {
	uint64_t ccid2 = cdbg.color_class_id(kmer2);
	
	if (ccid2 && ccid2 != ccid1) {
		if (edge(ccid1, ccid2, ccg).second)
			return;

		uint64_t num_samples = cdbg.get_num_samples();
		const colored_dbg::color_class cc1 = cdbg.get_color_class(ccid1);
		const colored_dbg::color_class cc2 = cdbg.get_color_class(ccid2);
		
		uint64_t hamming_distance = 0;
		for (uint64_t word_index = 0;
				 word_index < (num_samples + 63) / 64;
				 word_index ++) {
			uint64_t word1 = cc1.get_uint64(64 * word_index);
			uint64_t word2 = cc2.get_uint64(64 * word_index);
			hamming_distance += popcnt(word1 ^ word2);
		}

		add_edge(ccid1, ccid2, hamming_distance, ccg);
	}
}

static void add_zero_edge(colored_dbg &cdbg,
													ccGraph &ccg,
													uint64_t ccid) {
	uint64_t num_samples = cdbg.get_num_samples();
	const colored_dbg::color_class cc = cdbg.get_color_class(ccid);
	
	uint64_t hamming_distance = 0;
	for (uint64_t word_index = 0; word_index < (num_samples + 63) / 64; word_index ++) {
		uint64_t word = cc.get_uint64(word_index);
		hamming_distance += popcnt(word);
	}

	add_edge(0, ccid, hamming_distance, ccg);
}


int build_dbgccmst(DBGCCMSTOpts& opt)
{
  std::string prefix = opt.prefix;

  // Make sure the prefix is a full folder
  if (prefix.back() != '/') {
    prefix.push_back('/');
  }

  spdlog::logger* console = opt.console.get();
	console->info("Reading colored dbg from disk.");

	std::string dbg_file(prefix + mantis::CQF_FILE);
	std::string sample_file(prefix + mantis::SAMPLEID_FILE);
	std::vector<std::string> eqclass_files = mantis::fs::GetFilesExt(prefix.c_str(),
                                                                   mantis::EQCLASS_FILE);
	colored_dbg cdbg(dbg_file,
									 eqclass_files,
									 sample_file);

	const CQF<KeyObject> * cqf = cdbg.get_cqf();
	uint64_t num_samples = cdbg.get_num_samples();
	uint64_t num_ccs = cdbg.get_num_bitvectors();

	ccGraph ccg(num_ccs+1);
	for (auto it = cqf->begin(0); !it.done(); ++it) {
		auto ko = *it;
		uint64_t kmer1 = HashUtil::hash_64i(ko.key, BITMASK(cqf->keybits()));
		uint64_t ccid1 = ko.count;

		for (uint64_t b : {0, 1, 2, 3}) {
			uint64_t kmer_right = (kmer1 << 2 | b) & BITMASK(cqf->keybits());
			uint64_t kmer_left = (kmer1 >> 2) | (b << (cqf->keybits()-2));
			add_edge(cdbg, ccg, ccid1, kmer_right);
			add_edge(cdbg, ccg, ccid1, kmer_left);
		}
	}

	for (uint64_t ccid = 1; ccid < num_ccs + 1; ccid++)
		add_zero_edge(cdbg, ccg, ccid);


	printf("Graph:\n");
	auto pit = edges(ccg);
	for (auto it = pit.first; it != pit.second; ++it)
		printf("%ld %ld %ld\n",
					 source(*it, ccg),
					 target(*it, ccg),
					 get(edge_weight_t(), ccg, *it));


	std::vector<graph_traits < ccGraph >::edge_descriptor> mst;
	kruskal_minimum_spanning_tree(ccg, std::back_inserter(mst));
	
	printf("MST:\n");
	for (const auto & edg : mst)
		printf("%ld %ld %ld\n",
					 source(edg, ccg),
					 target(edg, ccg),
					 get(edge_weight_t(), ccg, edg));

	uint64_t total_weight = 0;
	for (const auto & edg : mst)
		total_weight += get(edge_weight_t(), ccg, edg);
	printf("Total MST weight: %ld\n",
				 total_weight);
	
	return EXIT_SUCCESS;
}
