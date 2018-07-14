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

using namespace boost;
typedef adjacency_list < vecS, vecS, undirectedS, no_property, property < edge_weight_t, uint64_t > > ccGraph;

// Check whether kmer2 exists and, if it does, add the corresponding
// edge to the color class graph.
//
// Assumes: kmer1 exists with color id ccid1.
static void add_edge(const CQF<KeyObject> *cqf,
										 uint64_t num_samples,
										 const BitVectorRRR & eqclasses,
										 ccGraph &ccg,
										 uint64_t ccid1,
										 uint64_t kmer2) {
	uint64_t kmer2rev = Kmer::reverse_complement(kmer2);
	if (Kmer::compare_kmers(kmer2rev, kmer2))
		kmer2 = kmer2rev;
		
	KeyObject qo(HashUtil::hash_64(kmer2, BITMASK(2*K)), 0, 0);
	uint64_t ccid2 = cqf->query(qo);
	if (ccid2 && ccid2 != ccid1) {
		if (edge(ccid1, ccid2, ccg).second)
			return;
		
		uint64_t hamming_distance = 0;
		for (uint64_t word_index = 0; word_index < (num_samples + 63) / 64; word_index ++) {
			uint64_t len = std::min(64UL, num_samples - 64 * word_index);
			uint64_t word1 = eqclasses.get_int((ccid1 - 1) * num_samples, len);
			uint64_t word2 = eqclasses.get_int((ccid2 - 1) * num_samples, len);
			hamming_distance += popcnt(word1 ^ word2);
		}

		add_edge(ccid1, ccid2, hamming_distance, ccg);
	}
}

static void add_zero_edge(uint64_t num_samples,
													const BitVectorRRR & eqclasses,
													ccGraph &ccg,
													uint64_t ccid) {
	uint64_t hamming_distance = 0;
	for (uint64_t word_index = 0; word_index < (num_samples + 63) / 64; word_index ++) {
		uint64_t len = std::min(64UL, num_samples - 64 * word_index);
		uint64_t word = eqclasses.get_int((ccid - 1) * num_samples, len);
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

	std::string cqf_file(prefix + CQF_FILE);
	std::string eqclass_file(prefix + EQCLASS_FILE);
	std::string sample_file(prefix + SAMPLEID_FILE);
	ColoredDbg<SampleObject<CQF<KeyObject>*>, KeyObject> cdbg(cqf_file,
																														eqclass_file,
																														sample_file);
  console->info("Read colored dbg with {} k-mers and {} color classes",
                cdbg.get_cqf()->size(),
								cdbg.get_bitvector().bit_size() / cdbg.get_num_samples());


	const CQF<KeyObject> *cqf = cdbg.get_cqf();
	const BitVectorRRR ccTable = cdbg.get_bitvector();
	uint64_t num_samples = cdbg.get_num_samples();
	uint64_t num_ccs = ccTable.bit_size() / num_samples;

	ccGraph ccg(num_ccs+1);
	for (auto it = cqf->begin(0); !it.done(); ++it) {
		auto ko = *it;
		uint64_t kmer1 = HashUtil::hash_64i(ko.key, BITMASK(2*K));
		uint64_t ccid1 = ko.count;

		for (uint64_t b : {0, 1, 2, 3}) {
			uint64_t kmer_right = (kmer1 << 2 | b) & BITMASK(2*K);
			uint64_t kmer_left = (kmer1 >> 2) | (b << (2*K-2));
			add_edge(cqf, num_samples, ccTable, ccg, ccid1, kmer_right);
			add_edge(cqf, num_samples, ccTable, ccg, ccid1, kmer_left);
		}
	}

	for (uint64_t ccid = 1; ccid < num_ccs + 1; ccid++)
		add_zero_edge(num_samples, ccTable, ccg, ccid);


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
