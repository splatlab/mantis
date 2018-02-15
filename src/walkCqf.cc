/*
 * ===========================================================================
 *
 *       Filename:  kmer_query.cc
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  04/27/2016 08:48:26 AM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Prashant Pandey (), ppandey@cs.stonybrook.edu
 *   Organization:  Stony Brook University
 *
 * ===========================================================================
 */
#include <iostream>
#include <algorithm>
#include <cstring>
#include <vector>
#include <set>
#include <bitset>
#include <cassert>
#include <fstream>

#include <boost/thread/thread.hpp>
#include <boost/lockfree/queue.hpp>
#include <boost/lockfree/spsc_queue.hpp>
#include <boost/atomic.hpp>

#include <time.h>
#include <stdio.h>
#include <fcntl.h>
#include <unistd.h>
#include <sys/resource.h>
#include <sys/stat.h>
#include <sys/time.h>
#include <sys/mman.h>
#include <cmath>

#include "cqf.h"
#include "hashutil.h"
#include "kmer.h"


#define BITMASK(nbits) ((nbits) == 64 ? 0xffffffffffffffff : (1ULL << (nbits)) - 1ULL)

/* Print elapsed time using the start and end timeval */
void print_time_elapsed(string desc, struct timeval* start, struct timeval* end) 
{
	struct timeval elapsed;
	if (start->tv_usec > end->tv_usec) {
		end->tv_usec += 1000000;
		end->tv_sec--;
	}
	elapsed.tv_usec = end->tv_usec - start->tv_usec;
	elapsed.tv_sec = end->tv_sec - start->tv_sec;
	float time_elapsed = (elapsed.tv_sec * 1000000 + elapsed.tv_usec)/1000000.f;
	std::cout << desc << "Total Time Elapsed: " << to_string(time_elapsed) << " seconds" << std::endl;
}

void run_filter(std::string ds_file, 
				std::string out_file, 
				uint64_t cutoff, 
				uint64_t approximate_num_of_kmers_greater_than_cutoff) {
	struct timeval start, end;
	struct timezone tzp;
	
	//Initialize the QF
	gettimeofday(&start, &tzp);
	cout << "Reading the input CQF off disk" << std::endl;
	CQF<KeyObject> cqf(ds_file, false);
	cout << "Done loading cqf in time: ";
	gettimeofday(&end, &tzp);
	print_time_elapsed("", &start, &end);
	typename CQF<KeyObject>::Iterator it = cqf.begin(0);

	uint64_t quotientBits = std::ceil(std::log2(approximate_num_of_kmers_greater_than_cutoff))+1;
	std::cout << "quotientBits : " << quotientBits << "\n";
	std::cout << "keybits : " << cqf.keybits() << "\n";
	CQF<KeyObject> newCqf(quotientBits, cqf.keybits(), cqf.seed());
	gettimeofday(&start, &tzp);
	uint64_t cntr = 1;
	uint64_t insertedCntr = 0;
	while (!it.done()) {
		KeyObject k = *it;
		if (k.count >= cutoff) {
			k.count = 1;//cutoff;
			newCqf.insert(k);
			insertedCntr++;
		}
		++it;
		if (cntr % 1000000 == 0) {
			std::cout << cntr << " kmers passed, " 
					  << insertedCntr << " kmers inserted, " 
					  << newCqf.noccupied_slots() << " slots occupied, "
					  << " load factor: " 
					  << static_cast<double>(newCqf.noccupied_slots())/static_cast<double>(1ULL << quotientBits) << "\n"; 
		}
		cntr++;
	}
	newCqf.serialize(out_file);
	gettimeofday(&end, &tzp);
	print_time_elapsed("", &start, &end);

}

void run_list_kmers(std::string ds_file, 
				std::string out_file, 
				uint64_t cutoff) {
	struct timeval start, end;
	struct timezone tzp;
	
	//Initialize the QF
	gettimeofday(&start, &tzp);
	cout << "Reading the input CQF off disk" << std::endl;
	CQF<KeyObject> cqf(ds_file, false);
	cout << "Done loading cqf in time: ";
	gettimeofday(&end, &tzp);
	print_time_elapsed("", &start, &end);
	typename CQF<KeyObject>::Iterator it = cqf.begin(0);

	gettimeofday(&start, &tzp);
	std::ofstream fout(out_file, ios::out);
	while (!it.done()) {
		KeyObject k = *it;
		// kmers.push_back(HashUtil::hash_64i(k.key, BITMASK(cqf.keybits())));
		uint64_t kint = HashUtil::hash_64i(k.key, BITMASK(cqf.keybits()));
		//if (k.count >= cutoff) {
			std::string kstr = Kmer::int_to_str(kint);
			fout << kstr << "\t" << k.count << "\n";
		//}
		++it;
	}
	fout.close();
	gettimeofday(&end, &tzp);
	print_time_elapsed("", &start, &end);
}

/*
 * ===  FUNCTION  ============================================================
 *         Name:  main
 *  Description:  
 * ===========================================================================
 */

int main ( int argc, char *argv[] )
{
	std::cout << argc << " ....... \n";
	std::string command = argv[1];
	if (command != "filter" and command != "list_kmers") {
		std::cerr <<  "ERROR: command can only be filter or list_kmers.\n";
		exit(1);
	}

	std::string ds_file = argv[2];
	cout << ds_file << "\n";
	std::string out_file = argv[3];
	cout << out_file << "\n";
	uint64_t cutoff = stoi(argv[4]);
	uint64_t approximate_num_of_kmers_greater_than_cutoff = 0;
	if (command == "filter") { 
		if (argc < 6) {
			std::cerr << "ERROR: missing last argument for filter command\n";
			exit(1);
		}
		cout << argv[5] << "\n";
		approximate_num_of_kmers_greater_than_cutoff = stoull(argv[5]);
	}

	if (command == "filter") {
		run_filter(ds_file, out_file, cutoff, approximate_num_of_kmers_greater_than_cutoff);
	}
	else {
		run_list_kmers(ds_file, out_file, cutoff);
	}
	
	return EXIT_SUCCESS;
}				/* ----------  end of function main  ---------- */

