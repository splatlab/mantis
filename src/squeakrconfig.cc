/*
 * =====================================================================================
 *
 *         Author:  Prashant Pandey (), ppandey@cs.stonybrook.edu
 *   Organization:  Stony Brook University
 *
 * =====================================================================================
 */

#include <fstream>      // std::ifstream

#include "squeakrconfig.h"

namespace squeakr {

	int read_config(std::string file, squeakrconfig *config) {
		// seek to the end of the file and read the k-mer size
		std::ifstream squeakr_file(file, std::ofstream::in);
		uint64_t file_size = squeakr_file.tellg();
		squeakr_file.seekg(file_size - sizeof(squeakrconfig));
		squeakr_file.read((char*)config, sizeof(config));
		squeakr_file.close();
		if (config->endianness != ENDIANNESS) {
			return SQUEAKR_INVALID_ENDIAN;
		}
		if (config->version != INDEX_VERSION) {
			return SQUEAKR_INVALID_VERSION;
		}
		return 0;
	}
}
