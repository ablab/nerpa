cmake_minimum_required(VERSION 3.5)
project(NRPsMatcher)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++14 -O3")

include_directories(src)

add_subdirectory(src)
add_subdirectory(tests)

install(PROGRAMS "${CMAKE_CURRENT_SOURCE_DIR}/run_nrp_matcher.py"
        DESTINATION bin
        COMPONENT runtime)