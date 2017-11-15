

add_library(mantis_core
  hashutil.cc
  coloreddbg.cc
  query.cc
  util.cc
  bitvector.cc
  cqf/gqf.c)

target_include_directories(mantis_core PUBLIC
  $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/include>)

target_link_libraries(mantis_core
  sdsl
  pthread
  boost_system 
  boost_thread
  m
  z
  rt)

add_executable(mantis
  mantis.cc)

target_include_directories(mantis PUBLIC
  $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/include>)

target_link_libraries(mantis
  mantis_core)
