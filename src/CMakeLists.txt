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
        NRP/NRP.h
        NRP/NRPsMatch.cpp
        NRP/NRPBuilder.cpp
        NRP/NRPBuilder.h
        NRP/NRPCycle.h
        NRP/NRPCycle.cpp
        NRP/NRPLine.h
        NRP/NRPLine.cpp
        NRP/NRPtail.h
        NRP/NRPtail.cpp
        NRPGenerator/NRPGenerator.h
        NRPGenerator/NRPGenerator.cpp
        NRPGenerator/NRPGeneratorTriplet.h
        NRPGenerator/NRPGeneratorTriplet.cpp
        NormalizedMatch/NormalizedMatch.h
        NormalizedMatch/NormalizedMatch.cpp)


add_library(nrps_matcher_lib ${SOURCE_FILES})
add_executable(NRPsMatcher ${SOURCE_FILES})

install(TARGETS NRPsMatcher
        DESTINATION bin
        COMPONENT runtime)