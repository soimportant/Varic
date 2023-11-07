cmake_minimum_required(VERSION 3.10)

set(SUBMODULES_DIR ${CMAKE_CURRENT_LIST_DIR})

# spdlog
if(NOT TARGET spdlog)
  add_subdirectory(${CMAKE_CURRENT_LIST_DIR}/spdlog)
endif()

# edlib
set(EDLIB_ENABLE_INSTALL OFF CACHE BOOL "Generate the install target" FORCE)
set(EDLIB_BUILD_EXAMPLES OFF CACHE BOOL "Build examples" FORCE)
set(EDLIB_BUILD_UTILITIES OFF CACHE BOOL "Build utilities" FORCE)
add_subdirectory(${CMAKE_CURRENT_LIST_DIR}/edlib)


# gperftools
# refer https://github.com/openfheorg/openfhe-development/blob/main/CMakeLists.txt

include(ExternalProject)

set(GPERFTOOLS_DIR ${SUBMODULES_DIR}/gperftools)
ExternalProject_add(
  gperftools
  SOURCE_DIR ${GPERFTOOLS_DIR}
  CMAKE_ARGS -DCMAKE_INSTALL_PREFIX=${GPERFTOOLS_DIR}
  BUILD_ALWAYS 1
  INSTALL_COMMAND make install
)
add_library(tcmalloc SHARED IMPORTED GLOBAL)
set_target_properties(tcmalloc PROPERTIES IMPORTED_LOCATION ${GPERFTOOLS_DIR}/lib/libtcmalloc.so)

add_library(profiler SHARED IMPORTED GLOBAL)
set_target_properties(profiler PROPERTIES IMPORTED_LOCATION ${GPERFTOOLS_DIR}/lib/libprofiler.so)
