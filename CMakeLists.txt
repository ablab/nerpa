cmake_minimum_required(VERSION 3.5)
project(Nerpa)

# allow more human readable "if then else" constructs
set(CMAKE_ALLOW_LOOSE_LOOP_CONSTRUCTS TRUE)

# modules to find additional libraries (like zlib, boost and etc)
# Add path for custom modules
set(CMAKE_MODULE_PATH
        ${CMAKE_MODULE_PATH}
        "${CMAKE_CURRENT_SOURCE_DIR}/src/cmake"
        "${CMAKE_CURRENT_SOURCE_DIR}/src/cmake/Modules")

# Define various dirs
set(MAIN_SRC_DIR ${CMAKE_CURRENT_SOURCE_DIR})
set(BINARY_DIR ${CMAKE_CURRENT_BINARY_DIR})
set(TOOLS_BINARY_DIR ${BINARY_DIR}/bin)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${TOOLS_BINARY_DIR})

# Sanity check our source directory to make sure that we are not trying to
# generate an in-tree build
if ((CMAKE_SOURCE_DIR STREQUAL CMAKE_BINARY_DIR) OR
(CMAKE_SOURCE_DIR STREQUAL "${CMAKE_BINARY_DIR}/src"))
    message(FATAL_ERROR "In-source builds are not allowed.
Please create a directory and run cmake from there, passing the path
to this source directory as the last argument.
This process created the file `CMakeCache.txt' and the directory `CMakeFiles'.
Please delete them.")
endif()


# Handle options, if any
include(options)

# Handle commond deps
include(deps)

# Add flags if necessary
include(flags)

# Handle common includes
include(includes)

# Useful macros
include(macros)

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

install(DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/configs"
        DESTINATION share/nerpa
        COMPONENT runtime)

install(DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/test_dataset"
        DESTINATION share/nerpa
        COMPONENT runtime)

install(DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/external_tools"
        DESTINATION share/nerpa
        COMPONENT runtime)
