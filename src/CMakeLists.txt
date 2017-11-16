
# make a library libmantis_core that will contain
# most of the relevant API
add_library(mantis_core
  kmer.cc
  hashutil.cc
  query.cc
  util.cc
  bitvector.cc
  validatemantis.cc
  coloreddbg.cc
  MantisFS.cc
  cqf/gqf.c)

target_include_directories(mantis_core PUBLIC
  $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/include>)

# link libmantis_core with the required libraries
target_link_libraries(mantis_core
  sdsl
  pthread
  m
  z
  rt)

# The driver program is `mantis`, which is a front-end
# with the relevant sub-commands
add_executable(mantis
  mantis.cc)

target_include_directories(mantis PUBLIC
  $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/include>)

target_link_libraries(mantis
  mantis_core)
