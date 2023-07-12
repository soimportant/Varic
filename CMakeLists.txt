cmake_minimum_required(VERSION 3.22)
set(CMAKE_CXX_COMPILER /usr/bin/g++-12)
set(CMAKE_C_COMPILER /usr/bin/gcc-12)
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

project(thesis LANGUAGES CXX)
set(CMAKE_BUILD_TYPE DEBUG)

# include_directories(${CMAKE_CURRENT_LIST_DIR}/lib)
# include_directories(${CMAKE_CURRENT_LIST_DIR}/include)

# biovoltron
set(BIOVOLTRON_TESTS OFF)
add_subdirectory(${CMAKE_CURRENT_LIST_DIR}/include/biovoltron)
add_subdirectory(${CMAKE_CURRENT_LIST_DIR}/submodules)

find_package(Boost COMPONENTS program_options thread graph REQUIRED)
find_package(OpenMP REQUIRED)

link_libraries(spdlog spoa::spoa ${Boost_LIBRARIES} biovoltron edlib OpenMP::OpenMP_CXX)
add_compile_definitions(DATA_PATH="/mnt/ec/mammoth/yolkee/thesis/data"
                        RESULT_PATH="/mnt/ec/mammoth/yolkee/thesis/result")
add_compile_options(-Wall -Wextra -Wpedantic)
include_directories(${CMAKE_CURRENT_LIST_DIR}/submodules/edlib/edlib/include)

add_executable(main ${CMAKE_CURRENT_LIST_DIR}/src/main.cpp)
# add_executable(test ${CMAKE_CURRENT_LIST_DIR}/src/test.cpp)
add_executable(eval ${CMAKE_CURRENT_LIST_DIR}/src/eval.cpp)
target_compile_definitions(eval INTERFACE _GLIBCXX_DEBUG)
add_executable(eval2 ${CMAKE_CURRENT_LIST_DIR}/src/eval2.cpp)
target_compile_definitions(eval2 INTERFACE _GLIBCXX_DEBUG)
add_executable(test ${CMAKE_CURRENT_LIST_DIR}/src/test.cpp)
