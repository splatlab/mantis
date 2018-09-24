/*
 * ============================================================================
 *
 *        Authors:  Prashant Pandey <ppandey@cs.stonybrook.edu>
 *                  Rob Johnson <robj@vmware.com>   
 *                  Rob Patro (rob.patro@cs.stonybrook.edu)
 *
 * ============================================================================
 */

#ifndef _SQUEAKR_CONFIG_H_
#define _SQUEAKR_CONFIG_H_

#include <string>
#include <stdio.h>

namespace squeakr {
	constexpr const uint32_t SQUEAKR_INVALID_VERSION{1};
	constexpr const uint32_t SQUEAKR_INVALID_ENDIAN{2};
	constexpr const uint32_t INDEX_VERSION{2};
	constexpr const uint64_t ENDIANNESS{0x0102030405060708ULL};

	typedef struct __attribute__ ((__packed__)) squeakrconfig {
		uint64_t kmer_size;
		uint64_t cutoff;
		uint64_t contains_counts;
		uint64_t endianness;
		uint32_t version;
	} squeakrconfig;

	int read_config(std::string squeakr_file, squeakrconfig *config);
}
#endif
