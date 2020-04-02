//
// Created by Fatemeh Almodaresi on 2020-03-17.
//

#include <iostream>
#include <fstream>
#include <sstream>
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

constexpr static uint32_t classicMaxCnt = 128;


int recursiveConstruction(uint64_t cnt, uint64_t start, uint64_t end, std::vector<std::pair<uint64_t, std::string>>&cmds,
                          std::string& complete_input_list, std::string& execDir, std::string& inputDir, std::string& tmpOut,
                          uint64_t numThreads,
                          uint64_t level) {
    std::string sysCommand;
    std::string tmpIn = "squeakr" + std::to_string(start) + "_" + std::to_string(end);
    if (end-start <= classicMaxCnt) {
        sysCommand = "head -" + std::to_string(end) + " " + complete_input_list
                     + " | tail -" + std::to_string(end-start) + " > " + inputDir + tmpIn + ";";
        sysCommand += "/usr/bin/time " + execDir + "/mantis build -s 20 -i " + inputDir + tmpIn
                     + " -o " + tmpOut + tmpIn + " -t "
                     + std::to_string(numThreads) + " > " + tmpOut + std::to_string(level) + tmpIn + ".log 2>&1;";
        sysCommand += "ls -lh " + tmpOut + tmpIn + " >> " + tmpOut + std::to_string(level) + tmpIn + ".log";
        cmds.emplace_back(level, sysCommand);
        return EXIT_SUCCESS;
    }

    uint64_t mid = (start+end)/2;
    recursiveConstruction(cnt, start, mid, cmds, complete_input_list,
                          execDir, inputDir, tmpOut, numThreads, level+1);
    recursiveConstruction(cnt, mid, end, cmds, complete_input_list,
                          execDir, inputDir, tmpOut, numThreads, level+1);
    sysCommand = "/usr/bin/time " + execDir + "/mantis merge -t " + std::to_string(numThreads)
                 + " -i1 " + tmpOut  + "squeakr" + std::to_string(start) + "_" + std::to_string(mid)
                 + " -i2 " + tmpOut  + "squeakr" + std::to_string(mid) + "_" + std::to_string(end)
                 + " -o " + tmpOut  + tmpIn
                 + " > " + tmpOut + std::to_string(level) + tmpIn + ".log 2>&1;";
    sysCommand += "ls -lh " + tmpOut + tmpIn + " >> " + tmpOut + std::to_string(level) + tmpIn + ".log;";
    sysCommand += "rm -rf " + tmpOut + "squeakr" + std::to_string(start) + "_" + std::to_string(mid) + ";";
    sysCommand += "rm -rf " + tmpOut + "squeakr" + std::to_string(mid) + "_" + std::to_string(end);
    cmds.emplace_back(level, sysCommand);
}

int construct_mantis_by_merge_main(BuildOpts &opt) {
    spdlog::logger *console = opt.console.get();
    std::ifstream infile(opt.inlist);

    uint32_t num_samples{0};
    if (infile.is_open()) {
        std::string line;
        while (std::getline(infile, line)) { ++num_samples; }
        infile.clear();
        infile.close();
        console->info("Will build mantis index over {} input experiments.", num_samples);
    } else {
        console->error("Input file {} does not exist or could not be opened.", opt.inlist);
        std::exit(1);
    }

    // Allowing only a limited num of processes run at the same time, even if more threads are provided
    ctpl::thread_pool p(opt.numProcesses /* num of threads in the pool */);
    std::vector<std::future<void>> results;
    results.reserve(num_samples);
    uint32_t numThreads{opt.numthreads};
    std::string complete_input_list{opt.inlist};

    std::string execDir = mantis::fs::getExecutableDir();

    if(opt.out.back() != '/')	// Make sure it is a full directory.
        opt.out += '/';
    if (!mantis::fs::DirExists(opt.out.c_str())) {
        mantis::fs::MakeDir(opt.out.c_str());
    } else if (not mantis::fs::IsDirEmpty(opt.out.c_str())){
        std::string cmd = "rm -rf " + opt.out + "*";
        system(cmd.c_str());
    }

    auto tmpOut = opt.out + mantis::TEMP_DIR;
    if (!mantis::fs::DirExists(tmpOut.c_str())) {
        mantis::fs::MakeDir(tmpOut.c_str());
    }
    auto inputDir = tmpOut + "inputs/";
    if (!mantis::fs::DirExists(inputDir.c_str())) {
        mantis::fs::MakeDir(inputDir.c_str());
    }
    std::vector<std::pair<uint64_t, std::string>> cmds;
    cmds.reserve(num_samples);
    recursiveConstruction(num_samples, 0, num_samples, cmds, complete_input_list,
            execDir, inputDir, tmpOut, numThreads, 0);
    if (cmds.empty()) {
        std::cerr << "ERROR: No command to run.\n";
        std::exit(2);
    }

    std::sort(cmds.begin(), cmds.end(), [](auto &c1, auto &c2){
        return c1.first > c2.first;
    });
    uint64_t level = cmds.begin()->first;
    for (auto &level_cmd:cmds) {
        if (level != level_cmd.first) {
            for (auto &r : results) {
                r.get();
            }
            results.clear();
            console->info("Done with Level {}", level);
            level = level_cmd.first;
        }
        results.push_back(p.push([level_cmd](int) {
            std::stringstream ss(std::to_string(level_cmd.first) + ": " + level_cmd.second + "\n\n");
            std::cerr << ss.str();
            system(level_cmd.second.c_str());
        }));
    }
    for (auto &r : results) {
        r.get();
    }
    std::string command = "mv "
                          + opt.out + mantis::TEMP_DIR + "squeakr" + std::to_string(0) + "_" + std::to_string(num_samples) + " "
                          + opt.out + "squeakr" + std::to_string(0) + "_" + std::to_string(num_samples) + ";"
                          + "rm -r " + inputDir + ";"
                          + "mv " + opt.out + mantis::TEMP_DIR + " " + opt.out + "logs;";
    system(command.c_str());
    std::cerr << "\nALL DONE!\n";
}
