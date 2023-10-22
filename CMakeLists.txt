cmake_minimum_required(VERSION 3.22)
set(CMAKE_CXX_COMPILER /usr/bin/g++-12)
set(CMAKE_C_COMPILER /usr/bin/gcc-12)
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

project(thesis LANGUAGES CXX)

# include_directories(${CMAKE_CURRENT_LIST_DIR}/lib)
include_directories(${CMAKE_CURRENT_LIST_DIR}/include)

# biovoltron
set(BIOVOLTRON_TESTS OFF)
add_subdirectory(${CMAKE_CURRENT_LIST_DIR}/include/biovoltron)

# submodules
add_subdirectory(${CMAKE_CURRENT_LIST_DIR}/submodules)

# spoa
set(spoa_optimize_for_native ON CACHE BOOL "optimize spoa with -march=native")
add_subdirectory(${CMAKE_CURRENT_LIST_DIR}/include/spoa)

find_package(Boost COMPONENTS program_options thread graph system filesystem REQUIRED)
find_package(OpenMP REQUIRED)

link_libraries(spdlog fmt spoa::spoa ${Boost_LIBRARIES} biovoltron edlib OpenMP::OpenMP_CXX)
include_directories(${CMAKE_CURRENT_LIST_DIR}/submodules/edlib/edlib/include)

if (CMAKE_BUILD_TYPE STREQUAL "Debug")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O2 -g -Wall -Wextra -march=native")
  # use flag below for debugging about memory issue
  # set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O2 -g -Wall -Wextra -lcurl -fsanitize=undefined -fsanitize=address")
  add_compile_definitions(DEBUG)
else()
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Ofast -march=native")
endif()

message(STATUS "CMAKE_CXX_FLAGS: ${CMAKE_CXX_FLAGS}")
add_compile_definitions(DATA_PATH="/mnt/ec/mammoth/yolkee/thesis/data"
RESULT_PATH="/mnt/ec/mammoth/yolkee/thesis/result")

add_executable(main ${CMAKE_CURRENT_LIST_DIR}/src/main.cpp)
# add_executable(eval ${CMAKE_CURRENT_LIST_DIR}/src/eval.cpp)
# add_executable(eval2 ${CMAKE_CURRENT_LIST_DIR}/src/eval2.cpp)
# add_executable(test ${CMAKE_CURRENT_LIST_DIR}/src/test.cpp)
# add_executable(maf2fa ${CMAKE_CURRENT_LIST_DIR}/src/maf2fa.cpp)


