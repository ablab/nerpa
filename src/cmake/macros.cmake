# -*- cmake -*-

# based on http://stackoverflow.com/questions/7787823/cmake-how-to-get-the-name-of-all-subdirectories-of-a-directory
#
# usage: subdirs_list(SUBDIRS ${MY_CURRENT_DIR})
macro(subdirs_list result curdir)
    file(GLOB children RELATIVE ${curdir} ${curdir}/*)
    set(dirlist "")
    foreach(child ${children})
        if(IS_DIRECTORY ${curdir}/${child})
            list(APPEND dirlist ${child})
        endif()
    endforeach()
    set(${result} ${dirlist})
endmacro(subdirs_list)

# based on http://cpansearch.perl.org/src/CBUREL/Alien-SmokeQt-4.6.0.4/cmake/modules/MacroOptionalAddSubdirectory.cmake
# Copyright (c) 2007, Alexander Neundorf, <neundorf@kde.org>
#
# Redistribution and use is allowed according to the terms of the BSD license.
macro(add_subdirectory_if_exists _dir)
    get_filename_component(_fullPath ${_dir} ABSOLUTE)
    if(EXISTS ${_fullPath}/CMakeLists.txt)
        add_subdirectory(${_dir})
    endif()
endmacro(add_subdirectory_if_exists)