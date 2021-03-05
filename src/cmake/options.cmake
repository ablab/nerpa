# -*- cmake -*-

# Build type
set(DEFAULT_BUILD_TYPE "RelWithDebInfo" CACHE STRING "Default build type")
if (NOT CMAKE_BUILD_TYPE)
  message("Setting default build configuration: ${DEFAULT_BUILD_TYPE}")
  set(CMAKE_BUILD_TYPE "${DEFAULT_BUILD_TYPE}" CACHE STRING
      "Choose the type of build, options are: None Debug Release RelWithAsserts RelWithDebInfo."
      FORCE)
endif()

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
    # -static-libgcc is gcc option only
    if (NOT "${CMAKE_CXX_COMPILER_ID}" MATCHES "Clang")
      set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -static-libgcc")
    endif()
  else()
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -static")
    add_definitions(-static)
  endif()
endif()
