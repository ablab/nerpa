cmake_minimum_required(VERSION 3.5)
project(Nerpa)

# Define option for static / dynamic build.
option(NERPA_STATIC_BUILD "Link Nerpa statically" OFF)
if (NERPA_STATIC_BUILD)
    # it'll make cmake to find libraries archives, not dynamic link
    set(CMAKE_FIND_LIBRARY_SUFFIXES .a)
    set(LINK_SEARCH_START_STATIC TRUE)
    set(LINK_SEARCH_END_STATIC TRUE)
    # This is dirty hack to get rid of -Wl,-Bdynamic
    set(CMAKE_EXE_LINK_DYNAMIC_C_FLAGS "-Wl,-Bstatic")
    set(CMAKE_EXE_LINK_DYNAMIC_CXX_FLAGS "-Wl,-Bstatic")

    if (APPLE)
        set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -static-libgcc")
    else()
        set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -static")
        add_definitions(-static)
    endif()

    set(Boost_USE_STATIC_LIBS        ON)
    set(Boost_USE_STATIC_RUNTIME     ON)
endif()


set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++14 -O3 -fopenmp")

include_directories(src)

add_subdirectory(src)
add_subdirectory(tests)

install(PROGRAMS "${CMAKE_CURRENT_SOURCE_DIR}/nerpa.py"
        DESTINATION bin
        COMPONENT runtime)

install(PROGRAMS "${CMAKE_CURRENT_SOURCE_DIR}/nerpa_init.py"
        DESTINATION bin
        COMPONENT runtime)

install(DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/src/nerpa_pipeline"
        DESTINATION share/nerpa
        COMPONENT runtime)

install(FILES
        "${CMAKE_CURRENT_SOURCE_DIR}/resources/test_dataset/library.info.streptomedb"
        "${CMAKE_CURRENT_SOURCE_DIR}/resources/test_dataset/prediction.info"
        "${CMAKE_CURRENT_SOURCE_DIR}/resources/aminoacids.tsv"
        "${CMAKE_CURRENT_SOURCE_DIR}/resources/modifications.tsv"
        "${CMAKE_CURRENT_SOURCE_DIR}/resources/monomers.tsv"
        "${CMAKE_CURRENT_SOURCE_DIR}/resources/monomersLogP.tsv"
        "${CMAKE_CURRENT_SOURCE_DIR}/resources/AAmod.tsv"
        DESTINATION share/nerpa
        COMPONENT runtime)

install(DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/resources/test_dataset/predictions"
        DESTINATION share/nerpa
        COMPONENT runtime)

install(DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/resources/test_dataset/streptomedb"
        DESTINATION share/nerpa
        COMPONENT runtime)


