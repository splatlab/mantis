#
# Help from :
# https://rix0r.nl/blog/2015/08/13/cmake-guide/
#

cmake_minimum_required(VERSION 3.0 FATAL_ERROR)
project(mantis VERSION 0.1 LANGUAGES C CXX)

if (NH)
   set(ARCH "")
   message(STATUS "Compiling mantis without Haswell instructions")
else()
    set(ARCH "-msse4.2 -D__SSE4_2_")
endif()

set(CMAKE_C_FLAGS "-std=gnu11 -Wno-unused-result -Wno-strict-aliasing -Wno-unused-function -Wno-sign-compare -Wno-implicit-function-declaration")
set(CMAKE_CXX_FLAGS "-std=c++11 ${ARCH} -Wno-unused-result -Wno-strict-aliasing -Wno-unused-function -Wno-sign-compare ${CMAKE_CXX_FLAGS}")
message("CMAKE_CXX_FLAGS = ${CMAKE_CXX_FLAGS}")
#set(CMAKE_EXE_LINKER_FLAGS "-L/home/rob/cosmo/3rd_party_inst/lib/")

include(GNUInstallDirs)

add_subdirectory(src)

