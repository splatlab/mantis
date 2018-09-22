/*
 * ============================================================================
 *
 *         Author:  Prashant Pandey (), ppandey@cs.stonybrook.edu
 *   Organization:  Stony Brook University
 *
 * ============================================================================
 */

#ifndef _UTIL_H_
#define _UTIL_H_

#include <iostream>
#include <cstring>
#include <vector>
#include <cassert>
#include <fstream>

#include <inttypes.h>

#ifdef DEBUG
#define PRINT_DEBUG 1
#else
#define PRINT_DEBUG 0
#endif

#define DEBUG(x) do { \
	if (PRINT_DEBUG) { std::cerr << x << std::endl; } \
} while (0)

#define ERROR(x) do { \
	{ std::cerr << x << std::endl; } \
} while (0)

#define PRINT(x) do { \
	{ std::cout << x << std::endl; } \
} while (0)

#define DEBUG_CDBG(x) do { \
	  if (PRINT_DEBUG) { std::cerr << x << std::endl; } \
} while (0)

#define PRINT_CDBG(x) do { \
	  { std::cout << x << std::endl; } \
} while (0)

class LightweightLock {
	public:
		LightweightLock() { locked = 0; }

		/**
		 * Try to acquire a lock once and return even if the lock is busy.
		 * If spin flag is set, then spin until the lock is available.
		 */
		void lock()
		{
			while (__sync_lock_test_and_set(&locked, 1))
				while (locked);
		}

		void unlock(void)
		{
			__sync_lock_release(&locked);
			return;
		}

	private:
		volatile int locked;
};

std::string last_part(std::string str, char c);
std::string first_part(std::string str, char c);
/* Print elapsed time using the start and end timeval */
void print_time_elapsed(std::string desc, struct timeval* start, struct
												timeval* end);
#endif
