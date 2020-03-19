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

constexpr static uint32_t classicMaxCnt = 25;

int construct_mantis_by_merge_main(BuildOpts &opt) {
    spdlog::logger *console = opt.console.get();
    std::ifstream infile(opt.inlist);

    uint32_t num_samples{0};
    if (infile.is_open()) {
        std::string line;
        while (std::getline(infile, line)) { ++num_samples; }
        infile.clear();
        infile.close();
//        infile.seekg(0, std::ios::beg);
        console->info("Will build mantis index over {} input experiments.", num_samples);
    } else {
        console->error("Input file {} does not exist or could not be opened.", opt.inlist);
        std::exit(1);
    }

    ctpl::thread_pool p(opt.numProcesses /* num of threads in the pool */);
    int arr[5] = {0};
    std::vector<std::future<void>> results;
    results.reserve(opt.numProcesses*10);
//    std::future<void> result;
    uint32_t scntr{0}, cntr{0};
    std::string sysCommand, createListCommand, input_lists;

    std::string execDir = mantis::fs::getExecutableDir();

    if(opt.out.back() != '/')	// Make sure it is a full directory.
        opt.out += '/';

    auto tmpOut = opt.out + mantis::TEMP_DIR;
    if (!mantis::fs::DirExists(tmpOut.c_str())) {
        mantis::fs::MakeDir(tmpOut.c_str());
    }
    auto inputDir = tmpOut + "inputs/";
    if (!mantis::fs::DirExists(inputDir.c_str())) {
        mantis::fs::MakeDir(inputDir.c_str());
    }
    while (scntr < num_samples) {
        auto tailCntr = std::min(num_samples-scntr, classicMaxCnt);
        input_lists = "squeakr" + std::to_string(scntr) + "_" + std::to_string(scntr + tailCntr);
        createListCommand = "head -" + std::to_string(scntr+tailCntr) + " " + opt.inlist
                + " | tail -" + std::to_string(tailCntr) + " > " + inputDir + input_lists;
        sysCommand = "/usr/bin/time " + execDir + "/mantis build -s 20 -i " + inputDir + input_lists
                     + " -o " + tmpOut + input_lists + " -t "
                     + std::to_string(opt.numthreads) + " > " + tmpOut + input_lists + ".log 2>&1";
        results.push_back(p.push([createListCommand, sysCommand, cntr](int) {
            std::cerr << "\n" << cntr << "\n\t" << createListCommand << "\n\t" << sysCommand << "\n";
//            usleep(10000000);
            system(createListCommand.c_str());
            system(sysCommand.c_str());
        }));
        scntr += tailCntr;
        cntr++;
    }
    for (auto &r : results) {
        r.get();
    }
    std::cerr << "\n\n\n\nall done\n";
}