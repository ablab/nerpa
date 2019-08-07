cmake_minimum_required(VERSION 3.5)
project(NRPsMatcher)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++14 -O3")

include_directories(src)

add_subdirectory(src)
add_subdirectory(tests)

install(PROGRAMS "${CMAKE_CURRENT_SOURCE_DIR}/run_nrp_matcher.py"
        DESTINATION bin
        COMPONENT runtime)

install(FILES
        "${CMAKE_CURRENT_SOURCE_DIR}/resources/test_dataset/library.info.streptomedb"
        "${CMAKE_CURRENT_SOURCE_DIR}/resources/test_dataset/prediction.info"
        "${CMAKE_CURRENT_SOURCE_DIR}/resources/aminoacids.tsv"
        "${CMAKE_CURRENT_SOURCE_DIR}/resources/modifications.tsv"
        "${CMAKE_CURRENT_SOURCE_DIR}/resources/AAmod.tsv"
        DESTINATION share/nrpsmatcher
        COMPONENT runtime)

install(DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/resources/test_dataset/predictions"
        DESTINATION share/nerpa
        COMPONENT runtime)

install(DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/resources/test_dataset/streptomedb"
        DESTINATION share/nerpa
        COMPONENT runtime)


