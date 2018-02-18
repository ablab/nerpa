cmake_minimum_required(VERSION 3.8)
project(NRPsMatcher)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++14 -O3")

add_subdirectory(tests)
add_subdirectory(src)

install(PROGRAMS "${CMAKE_CURRENT_SOURCE_DIR}/run_nrp_matcher.py"
        DESTINATION bin
        COMPONENT runtime)