# -*- cmake -*-

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++14")
add_compile_options(-Wno-deprecated)

# Use libc++ with clang due to C++11 mode
if ("${CMAKE_CXX_COMPILER_ID}" MATCHES "Clang")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -stdlib=libc++")
  # Require libsupc++ on Linux
  if (UNIX AND NOT APPLE)
    set(SYSTEM_LIBRARIES "supc++")
  endif()
endif()

if (${CMAKE_BUILD_TYPE} STREQUAL "Debug")
  message("Making Debug Configuration...")

  add_compile_options(-g3)
  add_definitions(-D_GLIBCXX_DEBUG)
  set(SPADES_DEBUG_LOGGING)
else()
  message("Making Release Configuration...")

  if (${CMAKE_BUILD_TYPE} STREQUAL "RelWithDebInfo")
    add_compile_options(-g3)
  else()
    add_compile_options(-g0)
  endif()

  add_compile_options(-O2)
  #  add_compile_options(-march=native)
  if (${CMAKE_BUILD_TYPE} STREQUAL "RelWithAsserts" OR
      ${CMAKE_BUILD_TYPE} STREQUAL "RelWithDebInfo")
    add_definitions(-UNDEBUG)
  else()
    add_definitions(-DNDEBUG)
  endif()
endif()

# Handle OpenMP flags
if ("${CMAKE_CXX_COMPILER_ID}" MATCHES "Clang" OR NOT OPENMP_FOUND)
  add_compile_options(-Wno-unknown-pragmas)
else ()
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp")
endif()

# Setup warnings
add_compile_options(-Wall -Wextra -Wconversion -Wno-sign-conversion -Wno-long-long -Wwrite-strings)