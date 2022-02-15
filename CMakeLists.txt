cmake_minimum_required(VERSION 3.15)
project(IOE_SORW)

set(CMAKE_CXX_STANDARD 14)

add_executable(GraSorw main.cpp)
add_executable(Convert2Undir scripts/Dir2Undir.cpp)
add_executable(MetisPreProcess scripts/MetisPre.cpp)
add_executable(GenPPRStarts scripts/GenPPRStarts.cpp)
add_executable(Dorder scripts/Dorder.cpp)
add_executable(CrawlWebConvert scripts/CWConvert.cpp)

include_directories(src/includes)
include_directories(src/apps)

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE RelWithDebInfo)
endif()

set(CMAKE_CXX_FLAGS "-fopenmp -Wall -Wno-strict-aliasing")
set(CMAKE_CXX_FLAGS_DEBUG "-g")
set(CMAKE_CXX_FLAGS_RELEASE "-O3")
set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "-g -O3")