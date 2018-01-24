cmake_minimum_required(VERSION 3.3)

project(NRPsMatcher_tests)

include_directories(${gtest_SOURCE_DIR}/include ${gtest_SOURCE_DIR})

add_executable(runContainNRPsPartTest containNRPsPartTest.cpp)

target_link_libraries(runContainNRPsPartTest gtest gtest_main)
target_link_libraries(runContainNRPsPartTest nrps_matcher_lib)

add_subdirectory(lib/gtest-1.7.0)