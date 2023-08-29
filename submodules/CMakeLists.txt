cmake_minimum_required(VERSION 3.10)


# spdlog
if(NOT TARGET spdlog)
  add_subdirectory(${CMAKE_CURRENT_LIST_DIR}/spdlog)
endif()

# edlib
set(EDLIB_BUILD_EXAMPLES OFF)
add_subdirectory(${CMAKE_CURRENT_LIST_DIR}/edlib)