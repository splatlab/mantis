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
#include <cmath>

#include <time.h>
#include <stdio.h>
#include <fcntl.h>
#include <unistd.h>
#include <sys/resource.h>
#include <sys/stat.h>
#include <sys/time.h>
#include <sys/mman.h>

#include "BooPHF.h"
#include "cqf.h"
#include "hashutil.h"
#include "kmer.h"


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

/*
 * ===  FUNCTION  ============================================================
 *         Name:  main
 *  Description:  
 * ===========================================================================
 */
int main ( int argc, char *argv[] )
{
	std::string kmers_file = argv[1];
	struct timeval start, end;
	struct timezone tzp;
	uint64_t kmer;
	std::vector<uint64_t> kmers;

	srand(time(NULL));

	gettimeofday(&start, &tzp);
	cout << "Reading the kmers off disk" << std::endl;
	std::ifstream fin(kmers_file, ios::in | ios::binary);
	kmers.reserve(4000000000);
	while(!fin.eof()) {
		fin >> kmer;
		kmers.push_back(kmer);
	}
    fin.close();
	std::cout << "Done loading the kmers" << std::endl;
	gettimeofday(&end, &tzp);
	print_time_elapsed("", &start, &end);
	
	gettimeofday(&start, &tzp);
	std::cout << "Building the Boom!" << std::endl;
	typedef boomphf::SingleHashFunctor<uint64_t> hasher_t;
	typedef boomphf::mphf<uint64_t, hasher_t> boophf_t;

	auto nkeys = kmers.size();
	auto keyIt = boomphf::range(kmers.begin(), kmers.end());
	boophf_t* bphf = new boophf_t(nkeys, keyIt, 16, 3.5);
	std::cout << "mphf size = " << (bphf->totalBitSize() / 8) / std::pow(2, 20)
			            << "\n";
	gettimeofday(&end, &tzp);
	print_time_elapsed("", &start, &end);

	return EXIT_SUCCESS;
}				/* ----------  end of function main  ---------- */

