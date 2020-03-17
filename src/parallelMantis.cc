//
// Created by Fatemeh Almodaresi on 2020-03-17.
//

#include <iostream>
#include <fstream>
//#include <utility>
//#include <algorithm>
//#include <string>
//#include <vector>
//#include <map>
//#include <queue>
//#include <set>
//#include <unordered_set>
//#include <bitset>
//#include <cassert>
//#include <fstream>

//#include <time.h>
//#include <stdio.h>
//#include <stdlib.h>
//#include <stdlib.h>
//#include <fcntl.h>
//#include <unistd.h>
//#include <sys/resource.h>
//#include <sys/stat.h>
//#include <sys/time.h>
//#include <sys/mman.h>
//#include <openssl/rand.h>

//#include "sparsepp/spp.h"
//#include "tsl/sparse_map.h"


#include "MantisFS.h"
#include "ProgOpts.h"
#include "ctpl_stl.h"
//#include "coloreddbg.h"
//#include "squeakrconfig.h"
//#include "json.hpp"
//#include "mantis_utils.hpp"
//#include "mantisconfig.hpp"

int construct_mantis_by_merge_main(BuildOpts &opt) {
    spdlog::logger *console = opt.console.get();
    std::ifstream infile(opt.inlist);

    uint64_t num_samples{0};
    if (infile.is_open()) {
        std::string line;
        while (std::getline(infile, line)) { ++num_samples; }
        infile.clear();
        infile.seekg(0, std::ios::beg);
        console->info("Will build mantis index over {} input experiments.", num_samples);
    } else {
        console->error("Input file {} does not exist or could not be opened.", opt.inlist);
        std::exit(1);
    }

    ctpl::thread_pool p(2 /* two threads in the pool */);
    int arr[5] = {0};
    std::vector<std::future<void>> results(4);
    for (int i = 0; i < 8; ++i) { // for 8 iterations,
        for (int j = 0; j < 4; ++j) {
            results[j] = p.push([&arr, j](int){ arr[j] +=2; });
        }
        for (int j = 0; j < 4; ++j) {
            results[j].get();
        }
        auto pt = std::min_element(arr, arr + 4);
        arr[4] = *pt;
    }
}