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

list(APPEND CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake/Modules")

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

###
# Sanitizers BEGIN
###
set(ASAN_LIB "")
#if (ASAN_BUILD)
#   #-fsanitize=leak
#   list(APPEND PF_CPP_FLAGS "-fsanitize=address")
#   #list(APPEND PF_CPP_FLAGS "-fsanitize=undefined")
#   #set(CMAKE_LINK_FLAGS "-fsanitize=address")
#   #list(APPEND CMAKE_LINK_FLAGS "-fsanitize=undefined")
#   set(ASAN_LIB "asan")
#endif()
###
# Sanitizers END
###

find_package(OpenMP)

if(NOT TARGET OpenMP::OpenMP_CXX)
   find_package(Threads REQUIRED)
   add_library(OpenMP::OpenMP_CXX IMPORTED INTERFACE)
   set_property(TARGET OpenMP::OpenMP_CXX
           PROPERTY INTERFACE_COMPILE_OPTIONS ${OpenMP_CXX_FLAGS})
   # Only works if the same flag is passed to the linker; use CMake 3.9+ otherwise (Intel, AppleClang)
   set_property(TARGET OpenMP::OpenMP_CXX
           PROPERTY INTERFACE_LINK_LIBRARIES ${OpenMP_CXX_FLAGS} Threads::Threads)
endif()

if (NOT GAT_SOURCE_DIR)
   message("Top-level source directory variable not set externally; setting it to ${CMAKE_CURRENT_SOURCE_DIR}")
   set(GAT_SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR})
endif()


##
# Let us check the sha sum of our packages if we have the right tools
##
set(SHASUM ${CMAKE_CURRENT_SOURCE_DIR}/scripts/check_shasum.sh)

## Try and find TBB first
find_package(TBB 2018.0 COMPONENTS tbb tbbmalloc tbbmalloc_proxy)

##
#
# Fetch and build Intel's Threading Building Blocks library.
#
##
if((NOT TBB_FOUND) OR (TBB_FOUND AND (TBB_VERSION VERSION_LESS 2018.0)))

   set(TBB_WILL_RECONFIGURE TRUE)
   # Set the appropriate compiler
   if(CLANG)
      set(TBB_COMPILER "clang")
   else()
      set(TBB_COMPILER "gcc")
   endif()

   message("Build system will fetch and build Intel Threading Building Blocks")
   message("==================================================================")
   # These are useful for the custom install step we'll do later
   set(TBB_SOURCE_DIR ${GAT_SOURCE_DIR}/external/oneTBB-2020.1)
   set(TBB_INSTALL_DIR ${GAT_SOURCE_DIR}/external/install)

   if("${TBB_COMPILER}" STREQUAL "gcc")
      ## Don't know why it's a problem yet, but if we're using
      ## GCC, get rid of the DO_ITT_NOTIFY flag
      set(TBB_CXXFLAGS "${TBB_CXXFLAGS} -UDO_ITT_NOTIFY")
   endif()

   set(TBB_CXXFLAGS "${TBB_CXXFLAGS} ${CXXSTDFLAG}")

   ExternalProject_Add(libtbb
           DOWNLOAD_DIR ${CMAKE_CURRENT_SOURCE_DIR}/external
           DOWNLOAD_COMMAND curl -k -L https://github.com/oneapi-src/oneTBB/archive/v2020.1.tar.gz -o v2020.1.tar.gz &&
           ${SHASUM} 7c96a150ed22bc3c6628bc3fef9ed475c00887b26d37bca61518d76a56510971 v2020.1.tar.gz &&
           tar -xzvf v2020.1.tar.gz
           ##
           #URL https://github.com/01org/tbb/archive/2018_U3.tar.gz
           #DOWNLOAD_NAME 2018_U3.tar.gz
           #URL_HASH SHA1=d6cf16a42ece60aad6a722b369e1a2aa753347b4
           #TLS_VERIFY FALSE
           ##
           SOURCE_DIR ${TBB_SOURCE_DIR}
           INSTALL_DIR ${TBB_INSTALL_DIR}
           PATCH_COMMAND "${TBB_PATCH_STEP}"
           CONFIGURE_COMMAND ""
           BUILD_COMMAND make ${QUIET_MAKE} CXXFLAGS=${TBB_CXXFLAGS} lambdas=1 compiler=${TBB_COMPILER} cfg=release tbb_build_prefix=LIBS
           INSTALL_COMMAND sh -c "mkdir -p ${TBB_INSTALL_DIR}/include && mkdir -p ${TBB_INSTALL_DIR}/lib && cp ${TBB_SOURCE_DIR}/build/LIBS_release/*.${SHARED_LIB_EXTENSION}* ${TBB_INSTALL_DIR}/lib && cp -r ${TBB_SOURCE_DIR}/include/* ${TBB_INSTALL_DIR}/include"
           BUILD_IN_SOURCE 1
           )

   set(RECONFIG_FLAGS ${RECONFIG_FLAGS} -DTBB_WILL_RECONFIGURE=FALSE -DTBB_RECONFIGURE=TRUE)
   ExternalProject_Add_Step(libtbb reconfigure
           COMMAND ${CMAKE_COMMAND} ${CMAKE_CURRENT_SOURCE_DIR} ${RECONFIG_FLAGS}
           DEPENDEES install
           )

   set(FETCHED_TBB TRUE)
endif()

##
# If we're fetching tbb, we need to have dummy paths for these variables
# so that CMake won't complain
##
if(TBB_WILL_RECONFIGURE)
   set(TBB_INCLUDE_DIRS ${TBB_INSTALL_DIR}/include)
   set(TBB_INCLUDE_DIR ${TBB_INSTALL_DIR}/include)
   set(TBB_LIBRARY_DIRS ${TBB_INSTALL_DIR}/lib)
   set(TBB_LIB_DIR ${TBB_INSTALL_DIR}/lib)
   #set(TBB_LIBRARIES tbb tbbmalloc)
   set(TBB_LIBRARIES ${TBB_INSTALL_DIR}/lib/libtbb.${SHARED_LIB_EXTENSION}
           ${TBB_INSTALL_DIR}/lib/libtbbmalloc.${SHARED_LIB_EXTENSION}
           )
   message("TBB_INCLUDE_DIRS = ${TBB_INCLUDE_DIRS}")
   message("TBB_LIBRARY_DIRS = ${TBB_LIBRARY_DIRS}")
endif()

##
#  The libtbb reconfigure should force this code
#  to be run on the second configuration pass, where it should appropriately set the
#  TBB_INSTALL_DIR variable.
##
if(TBB_RECONFIGURE)
   unset(TBB_FOUND CACHE)
   unset(TBB_INSTALL_DIR CACHE)
   unset(CMAKE_PREFIX_PATH CACHE)
   set(CMAKE_PREFIX_PATH ${GAT_SOURCE_DIR}/external/install)
   set(TBB_INSTALL_DIR ${GAT_SOURCE_DIR}/external/install)
   message("TBB_INSTALL_DIR = ${TBB_INSTALL_DIR}")
   find_package(TBB 2018.0 COMPONENTS tbb tbbmalloc tbbmalloc_proxy)
endif()
message("TBB_LIBRARIES = ${TBB_LIBRARIES}")



set(MANTIS_C_WARN "-Wno-unused-result;-Wno-strict-aliasing;-Wno-unused-function;-Wno-sign-compare;-Wno-implicit-function-declaration")
set(MANTIS_CXX_WARN "-Wno-unused-result;-Wno-strict-aliasing;-Wno-unused-function;-Wno-sign-compare")
set(MANTIS_C_FLAGS "${ARCH_DEFS};${MANTIS_C_WARN};${PF_CPP_FLAGS};")
set(MANTIS_CXX_FLAGS "${ARCH_DEFS};${MANTIS_CXX_WARN};${PF_CPP_FLAGS};")

if (SDSL_INSTALL_PATH)
   message("Adding ${SDSL_INSTALL_PATH}/include to the include path")
   message("Adding ${SDSL_INSTALL_PATH}/build/lib to the build path")
   list(APPEND MANTIS_CXX_FLAGS "-I${SDSL_INSTALL_PATH}/include")
   list(APPEND MANTIS_CXX_FLAGS "-L${SDSL_INSTALL_PATH}/lib")
endif()

include_directories(include)
link_directories(lib)


include(GNUInstallDirs)

add_subdirectory(src)

