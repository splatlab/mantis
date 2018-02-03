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

/*
 * ===  FUNCTION  ============================================================
 *         Name:  main
 *  Description:  
 * ===========================================================================
 */
int main ( int argc, char *argv[] )
{
	std::string ds_file = argv[1];
	struct timeval start, end;
	struct timezone tzp;
	std::vector<uint64_t> kmers;

	srand(time(NULL));

	//Initialize the QF
	cout << "Reading the input CQF off disk" << std::endl;
	CQF<KeyObject> cqf(ds_file, false);
	typename CQF<KeyObject>::Iterator it = cqf.begin(0);

	gettimeofday(&start, &tzp);
	while (!it.done()) {
		KeyObject k = *it;
		kmers.push_back(HashUtil::hash_64i(k.key, BITMASK(cqf.keybits())));
		++it;
	}
	std::cout << "Done reading" << std::endl;
	gettimeofday(&end, &tzp);
	print_time_elapsed("", &start, &end);

	
	gettimeofday(&start, &tzp);
	std::ofstream fout("kmers.bin", ios::out | ios::binary);
	size_t num = kmers.size();
    size_t elemSize = sizeof(typename std::vector<uint64_t>::value_type);
    // We have to get rid of constness below, but this should be OK
    fout.write(reinterpret_cast<char*>(const_cast<uint64_t*>(kmers.data())),
                  num * elemSize);
    fout.close();
	gettimeofday(&end, &tzp);
	print_time_elapsed("", &start, &end);

	return EXIT_SUCCESS;
}				/* ----------  end of function main  ---------- */

