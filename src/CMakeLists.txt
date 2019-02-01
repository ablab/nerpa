include_directories(.)

add_subdirectory(Logger)
set(SOURCE_FILES main.cpp
        Aminoacid/Aminoacid.cpp
        Aminoacid/Aminoacid.h
        NRPsPrediction/AminoacidPrediction.cpp
        NRPsPrediction/AminoacidPrediction.h
        NRPsPrediction/NRPsPart.cpp
        NRPsPrediction/NRPsPart.h
        NRPsPrediction/NRPsPrediction.cpp
        NRPsPrediction/NRPsPrediction.h
        NRP/NRP.cpp
        NRP/NRP.h
        Matcher/NRPsMatch.cpp
        NRP/NRPBuilder.cpp
        NRP/NRPBuilder.h
        NRP/NRPCycle.h
        NRP/NRPCycle.cpp
        NRP/NRPLine.h
        NRP/NRPLine.cpp
        NRP/NRPtail.h
        NRP/NRPtail.cpp
        Matcher/Matcher.h
        Matcher/Matcher.cpp
        Matcher/Score/Score.cpp
        Matcher/Score/Score.h
        Matcher/Segment.h
        NRPsPrediction/Builders/PredictionBuilderBase.h
        NRPsPrediction/Builders/PredictionBuilderBase.cpp
        NRPsPrediction/Builders/Nrpspredictor2Builder.h
        NRPsPrediction/Builders/Nrpspredictor2Builder.cpp
        Aminoacid/Formula.h
        Aminoacid/Modification.h
        Aminoacid/Modification.cpp)

add_library(nrps_matcher_lib ${SOURCE_FILES})
add_executable(NRPsMatcher main.cpp)

target_link_libraries(NRPsMatcher logger_lib nrps_matcher_lib)

install(TARGETS NRPsMatcher
        DESTINATION bin
        COMPONENT runtime)