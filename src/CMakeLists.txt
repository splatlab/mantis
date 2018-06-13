
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
add_executable(monochromatic_component_iterator
				monochromatic_component_iterator.cc)
add_executable(walkEqcls walkEqcls.cc)

target_include_directories(mantis PUBLIC
  $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/include>)

target_include_directories(monochromatic_component_iterator PUBLIC
  $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/include>)

target_include_directories(walkEqcls PUBLIC
  $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/include> $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/include/SIMDCompressionAndIntersection>)

target_link_libraries(mantis
  mantis_core)
target_link_libraries(monochromatic_component_iterator
				mantis_core)
target_link_libraries(walkEqcls
				mantis_core)
