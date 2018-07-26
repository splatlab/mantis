
# make a library libmantis_core that will contain
# most of the relevant API
add_library(mantis_core STATIC
  kmer.cc
  hashutil.cc
  query.cc
  util.cc
  bitvector.cc
  validatemantis.cc
  coloreddbg.cc
  MantisFS.cc
  cqf/gqf.c)

set(MANTIS_DEBUG_CFLAGS "${MANTIS_C_FLAGS};-g")
set(MANTIS_DEBUG_CXXFLAGS "${MANTIS_CXX_FLAGS};-g")
set(MANTIS_RELEASE_CFLAGS "${MANTIS_C_FLAGS};-O3;-march=native")
set(MANTIS_RELEASE_CXXFLAGS "${MANTIS_CXX_FLAGS};-O3;-march=native")

target_include_directories(mantis_core PUBLIC $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/include>)
target_compile_options(mantis_core PUBLIC "$<$<AND:$<CONFIG:DEBUG>,$<COMPILE_LANGUAGE:C>>:${MANTIS_DEBUG_CFLAGS}>")
target_compile_options(mantis_core PUBLIC "$<$<AND:$<CONFIG:DEBUG>,$<COMPILE_LANGUAGE:CXX>>:${MANTIS_DEBUG_CXXFLAGS}>")
target_compile_options(mantis_core PUBLIC "$<$<AND:$<CONFIG:RELEASE>,$<COMPILE_LANGUAGE:C>>:${MANTIS_RELEASE_CFLAGS}>")
target_compile_options(mantis_core PUBLIC "$<$<AND:$<CONFIG:RELEASE>,$<COMPILE_LANGUAGE:CXX>>:${MANTIS_RELEASE_CXXFLAGS}>")
target_compile_definitions(mantis_core PUBLIC "${ARCH_DEFS}")
set_property(TARGET mantis_core PROPERTY INTERPROCEDURAL_OPTIMIZATION True)

# link libmantis_core with the required libraries
target_link_libraries(mantis_core
  Threads::Threads
  m
  z
  rt
  sdsl)

# The driver program is `mantis`, which is a front-end
# with the relevant sub-commands
add_executable(mantis mantis.cc)
target_include_directories(mantis PUBLIC $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/include>)
target_link_libraries(mantis mantis_core)
target_compile_options(mantis PUBLIC "$<$<AND:$<CONFIG:DEBUG>,$<COMPILE_LANGUAGE:C>>:${MANTIS_DEBUG_CFLAGS}>")
target_compile_options(mantis PUBLIC "$<$<AND:$<CONFIG:DEBUG>,$<COMPILE_LANGUAGE:CXX>>:${MANTIS_DEBUG_CXXFLAGS}>")
target_compile_options(mantis PUBLIC "$<$<AND:$<CONFIG:RELEASE>,$<COMPILE_LANGUAGE:C>>:${MANTIS_RELEASE_CFLAGS}>")
target_compile_options(mantis PUBLIC "$<$<AND:$<CONFIG:RELEASE>,$<COMPILE_LANGUAGE:CXX>>:${MANTIS_RELEASE_CXXFLAGS}>")
target_compile_definitions(mantis PUBLIC "${ARCH_DEFS}")
set_property(TARGET mantis PROPERTY INTERPROCEDURAL_OPTIMIZATION True)

add_executable(build_eq_graph build_eq_graph.cc)
target_include_directories(build_eq_graph PUBLIC $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/include>)
target_link_libraries(build_eq_graph mantis_core)
target_compile_options(build_eq_graph PUBLIC "$<$<AND:$<CONFIG:DEBUG>,$<COMPILE_LANGUAGE:C>>:${MANTIS_DEBUG_CFLAGS}>")
target_compile_options(build_eq_graph PUBLIC "$<$<AND:$<CONFIG:DEBUG>,$<COMPILE_LANGUAGE:CXX>>:${MANTIS_DEBUG_CXXFLAGS}>")
target_compile_options(build_eq_graph PUBLIC "$<$<AND:$<CONFIG:RELEASE>,$<COMPILE_LANGUAGE:C>>:${MANTIS_RELEASE_CFLAGS}>")
target_compile_options(build_eq_graph PUBLIC "$<$<AND:$<CONFIG:RELEASE>,$<COMPILE_LANGUAGE:CXX>>:${MANTIS_RELEASE_CXXFLAGS}>")
target_compile_definitions(build_eq_graph PUBLIC "${ARCH_DEFS}")
set_property(TARGET build_eq_graph PROPERTY INTERPROCEDURAL_OPTIMIZATION True)

# TODO: look more into why this is necessary
if (SDSL_INSTALL_PATH)
   set_property(TARGET mantis APPEND_STRING PROPERTY LINK_FLAGS "-L${SDSL_INSTALL_PATH}/lib")
   set_property(TARGET mantis_core APPEND_STRING PROPERTY LINK_FLAGS "-L${SDSL_INSTALL_PATH}/lib")
endif()

install(TARGETS mantis
        RUNTIME DESTINATION bin)

add_executable(monochromatic_component_iterator
		monochromatic_component_iterator.cc)
target_include_directories(monochromatic_component_iterator PUBLIC
  $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/include>)
target_link_libraries(monochromatic_component_iterator
				mantis_core)
set_property(TARGET monochromatic_component_iterator PROPERTY INTERPROCEDURAL_OPTIMIZATION True)


add_executable(msf MSF.cc)
target_include_directories(msf PUBLIC
				  $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/include>)
target_link_libraries(msf
				                mantis_core)
set_property(TARGET msf PROPERTY INTERPROCEDURAL_OPTIMIZATION True)

add_executable(loadColorCls validateMSF.cc)
target_include_directories(loadColorCls PUBLIC
		$<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/include>)
target_link_libraries(loadColorCls
		mantis_core)
set_property(TARGET loadColorCls PROPERTY INTERPROCEDURAL_OPTIMIZATION True)

add_executable(mstBoost MST_boost.cc)
		target_include_directories(mstBoost PUBLIC
						                  $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/include>)
target_link_libraries(mstBoost mantis_core)
set_property(TARGET mstBoost PROPERTY INTERPROCEDURAL_OPTIMIZATION True)

#add_executable(walkEqcls walkEqcls.cc)
#target_include_directories(walkEqcls PUBLIC
#		$<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/include> $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/include/SIMDCompressionAndIntersection>)
#find_library(compression_library NAMES SIMDCompressionAndIntersection HINTS "${CMAKE_SOURCE_DIR}/lib")
#target_link_libraries(walkEqcls
#				mantis_core ${compression_library})
