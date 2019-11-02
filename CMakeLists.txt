#
# Help from :
# https://rix0r.nl/blog/2015/08/13/cmake-guide/
#

cmake_minimum_required(VERSION 3.5 FATAL_ERROR)
include(ExternalProject)

project(mantis VERSION 0.2 LANGUAGES C CXX)
if (NOT CMAKE_BUILD_TYPE)
   set (CMAKE_BUILD_TYPE "Release")
endif()

# We require C++17
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_C_STANDARD 11)
set(CMAKE_C_STANDARD_REQUIRED ON)


set(THREADS_PREFER_PTHREAD_FLAG ON)
find_package(Threads REQUIRED)

if (NH)
   set(ARCH_FLAGS "")
   set(ARCH_DEFS "")
   message(STATUS "Compiling mantis without Haswell instructions")
else()
   set(ARCH_FLAGS "")
   set(ARCH_DEFS "-D__SSE4_2_")
endif()

set(MANTIS_C_WARN "-Wno-unused-result;-Wno-strict-aliasing;-Wno-unused-function;-Wno-sign-compare;-Wno-implicit-function-declaration")
set(MANTIS_CXX_WARN "-Wno-unused-result;-Wno-strict-aliasing;-Wno-unused-function;-Wno-sign-compare")
set(MANTIS_C_FLAGS "${ARCH_DEFS};${MANTIS_C_WARN}")
set(MANTIS_CXX_FLAGS "${ARCH_DEFS};${MANTIS_CXX_WARN}")

if (SDSL_INSTALL_PATH)
   message("Adding ${SDSL_INSTALL_PATH}/include to the include path")
   message("Adding ${SDSL_INSTALL_PATH}/build/lib to the build path")
   list(APPEND MANTIS_CXX_FLAGS "-I${SDSL_INSTALL_PATH}/include")
   list(APPEND MANTIS_CXX_FLAGS "-L${SDSL_INSTALL_PATH}/lib")
endif()
#[[
ExternalProject_Add(libbf
        GIT_REPOSITORY https://github.com/mavam/libbf.git
        GIT_TAG        master
        UPDATE_COMMAND ""
        UPDATE_DISCONNECTED 1
        BUILD_IN_SOURCE TRUE
        DOWNLOAD_DIR ${CMAKE_CURRENT_SOURCE_DIR}/external/libbf
        SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/external/libbf
        INSTALL_DIR ${CMAKE_CURRENT_SOURCE_DIR}/external/install
        CONFIGURE_COMMAND ./configure
        BUILD_COMMAND     make
#        CXXFLAGS='-std=c++14'
        INSTALL_COMMAND   mkdir -p <INSTALL_DIR>/lib && mkdir -p <INSTALL_DIR>/include &&
        cp build/lib/libbf.so <INSTALL_DIR>/lib &&
        cp -r bf <INSTALL_DIR>/include
        )
]]

#add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/external/libbf)
#set(LIB_FLAGS "-L ${CMAKE_CURRENT_SOURCE_DIR}/external/lib")
set(LIB_FLAGS "")
include_directories(include #[[${CMAKE_CURRENT_SOURCE_DIR}/external/install/include]])
link_directories(lib #[[${CMAKE_CURRENT_SOURCE_DIR}/external/install/lib]])


include(GNUInstallDirs)

add_subdirectory(src)

