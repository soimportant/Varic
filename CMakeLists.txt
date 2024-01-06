cmake_minimum_required(VERSION 3.22)

project(thesis LANGUAGES C CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# if not define CMAKE_BUILD_TYPE, set it to Release
if (NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE "Release")
endif()

message("Build type: ${CMAKE_BUILD_TYPE}")
set(CMAKE_CXX_FLAGS_RELEASE "-O2")
set(CMAKE_CXX_FLAGS_DEBUG "-O2 -g")
# set(CMAKE_CXX_FLAGS_DEBUG "-O2 -g -fsanitize=address -fsanitize=undefined")

# boost
# set(BOOST_ROOT "/mnt/ec/ness/yolkee/miniconda3/envs/bio/include/boost/")
# set(BOOST_LIBRARYDIR "/mnt/ec/ness/yolkee/miniconda3/envs/bio/include")
find_package(Boost COMPONENTS program_options thread graph system filesystem REQUIRED)
link_libraries(${Boost_LIBRARIES})
include_directories(${Boost_INCLUDE_DIRS})

# Openmp
find_package(OpenMP REQUIRED)
link_libraries(OpenMP::OpenMP_C OpenMP::OpenMP_CXX)

# biovoltron
set(BIOVOLTRON_TESTS OFF)
add_subdirectory(${CMAKE_CURRENT_LIST_DIR}/include/biovoltron)

# submodules
add_subdirectory(${CMAKE_CURRENT_LIST_DIR}/submodules)

# spoa
add_subdirectory(${CMAKE_CURRENT_LIST_DIR}/include/spoa)

link_libraries(spdlog spoa::spoa ${Boost_LIBRARIES} biovoltron edlib OpenMP::OpenMP_C OpenMP::OpenMP_CXX tcmalloc profiler)

include_directories(${CMAKE_CURRENT_LIST_DIR}/include)

add_compile_definitions(DATA_PATH="/mnt/ec/mammoth/yolkee/thesis/data"
RESULT_PATH="/mnt/ec/mammoth/yolkee/thesis/result"
TMP_PATH="/mnt/ec/ness/yolkee/thesis/tmp")

set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -Wall -Wextra")

if ("${CMAKE_CXX_COMPILER_ID}" MATCHES "Clang")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fuse-ld=lld")
  # raise compile error on cannot find header, figure it out
  # set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -stdlib=libc++")
endif()

add_executable(main ${CMAKE_CURRENT_LIST_DIR}/src/main.cpp)
# add_executable(eval ${CMAKE_CURRENT_LIST_DIR}/src/eval.cpp)
# add_executable(eval2 ${CMAKE_CURRENT_LIST_DIR}/src/eval2.cpp)
add_executable(test ${CMAKE_CURRENT_LIST_DIR}/src/test.cpp)
add_executable(maf2fa ${CMAKE_CURRENT_LIST_DIR}/src/maf2fa.cpp)


