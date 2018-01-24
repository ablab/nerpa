include_directories(.)

set(SOURCE_FILES main.cpp
        Aminoacid/Aminoacids.cpp
        Aminoacid/Aminoacids.h
        NRPsPrediction/AminoacidPrediction.cpp
        NRPsPrediction/AminoacidPrediction.h
        NRPsPrediction/NRPsPart.cpp
        NRPsPrediction/NRPsPart.h
        NRPsPrediction/NRPsPrediction.cpp
        NRPsPrediction/NRPsPrediction.h
        NRP/NRP.cpp
        NRP/NRP.h)

add_library(nrps_matcher_lib ${SOURCE_FILES})
add_executable(NRPsMatcher ${SOURCE_FILES})