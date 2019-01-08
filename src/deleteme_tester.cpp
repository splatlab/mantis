//
// Created by Fatemeh Almodaresi on 2019-01-08.
//
#include "gqf_cpp.h"

int main ( int argc, char *argv[] ) {
    std::string inputfile = argv[1];
    std::string cqf_file(inputfile);
    CQF<KeyObject> cqf(cqf_file, CQF_FREAD);
    cqf.free();
    int tmp;
    std::cerr << "\nhere\n";
    std::cin >> tmp;
    std::cerr << "\n" << tmp << "\n";
}