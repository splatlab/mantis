
# make a library libmantis_core that will contain
# most of the relevant API
add_library(mantis_core STATIC
		kmer.cc
		query.cc
		mstQuery.cc
        validateMST.cc
		util.cc
  		validatemantis.cc
  		coloreddbg.cc
		canonicalKmer.cc
  		mst.cc
		stat.cc
  		MantisFS.cc
  		squeakrconfig.cc
  		gqf/gqf.c
  		gqf/gqf_file.c
  		gqf/hashutil.c)

set(MANTIS_DEBUG_CFLAGS "${MANTIS_C_FLAGS};-g")
set(MANTIS_DEBUG_CXXFLAGS "${MANTIS_CXX_FLAGS};-g")
set(MANTIS_RELEASE_CFLAGS "${MANTIS_C_FLAGS};-O3")
set(MANTIS_RELEASE_CXXFLAGS "${MANTIS_CXX_FLAGS};-O3")

target_include_directories(mantis_core PUBLIC $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/include>)
target_compile_options(mantis_core PUBLIC "$<$<AND:$<CONFIG:DEBUG>,$<COMPILE_LANGUAGE:C>>:${MANTIS_DEBUG_CFLAGS}>")
target_compile_options(mantis_core PUBLIC "$<$<AND:$<CONFIG:DEBUG>,$<COMPILE_LANGUAGE:CXX>>:${MANTIS_DEBUG_CXXFLAGS}>")
target_compile_options(mantis_core PUBLIC "$<$<AND:$<CONFIG:RELEASE>,$<COMPILE_LANGUAGE:C>>:${MANTIS_RELEASE_CFLAGS}>")
target_compile_options(mantis_core PUBLIC "$<$<AND:$<CONFIG:RELEASE>,$<COMPILE_LANGUAGE:CXX>>:${MANTIS_RELEASE_CXXFLAGS}>")
target_compile_definitions(mantis_core PUBLIC "${ARCH_DEFS}")

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

#add_executable(estimateNumOfKners estimateNumOfKmers.cc)
#target_include_directories(estimateNumOfKners PUBLIC $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/include>)
#target_link_libraries(estimateNumOfKners mantis_core)

# TODO: look more into why this is necessary
if (SDSL_INSTALL_PATH)
   set_property(TARGET mantis APPEND_STRING PROPERTY LINK_FLAGS "-L${SDSL_INSTALL_PATH}/lib")
   set_property(TARGET mantis_core APPEND_STRING PROPERTY LINK_FLAGS "-L${SDSL_INSTALL_PATH}/lib")
endif()

install(TARGETS mantis
        RUNTIME DESTINATION bin)
