# -*- cmake -*-

if ("${CMAKE_CXX_COMPILER_ID}" MATCHES "GNU")
    # Require at least gcc 4.9
    if (CMAKE_CXX_COMPILER_VERSION VERSION_LESS 4.9)
        message(FATAL_ERROR "gcc version 4.9 or later is required")
    endif()
elseif ("${CMAKE_CXX_COMPILER_ID}" MATCHES "Clang")
    if (CMAKE_CXX_COMPILER_VERSION VERSION_LESS 3.2)
        message(FATAL_ERROR "clang version 3.2 or later is required")
    endif()
else()
    message(WARNING "Unsupported compiler is detected. Compilation was not tested on it and may fail")
endif()

find_package(OpenMP)
