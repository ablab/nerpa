include_directories(.)

add_subdirectory(Logger)
set(SOURCE_FILES main.cpp
        utils/openmp_wrapper.h
        utils/cxxopts.hpp
        Aminoacid/Aminoacid.cpp
        Aminoacid/Aminoacid.h
        NRPsPrediction/AAdomainPrediction.cpp
        NRPsPrediction/AAdomainPrediction.h
        NRPsPrediction/OrfPrediction.cpp
        NRPsPrediction/OrfPrediction.h
        NRPsPrediction/BgcPrediction.cpp
        NRPsPrediction/BgcPrediction.h
        NRP/NRP.cpp
        NRP/NRP.h
        Matcher/NRPsMatch.cpp
        NRP/NRPBuilder.cpp
        NRP/NRPBuilder.h
        NRP/NRPCycle.h
        NRP/NRPCycle.cpp
        NRP/NRPLine.h
        NRP/NRPLine.cpp
        NRP/NRPBranch.h
        NRP/NRPBranch.cpp
        Matcher/OrderedGenesMatcher.h
        Matcher/OrderedGenesMatcher.cpp
        Matcher/MatcherBase.h
        Matcher/Score/Score.cpp
        Matcher/Score/Score.h
        Matcher/Segment.h
        NRPsPrediction/Builders/PredictionBuilderBase.h
        NRPsPrediction/Builders/PredictionBuilderBase.cpp
        NRPsPrediction/Builders/Nrpspredictor2Builder.h
        NRPsPrediction/Builders/Nrpspredictor2Builder.cpp
        Aminoacid/Formula.h
        Aminoacid/Formula.cpp
        Aminoacid/Modification.h
        Aminoacid/Modification.cpp
        Aminoacid/AminoacidInfo.h
        Aminoacid/AminoacidInfo.cpp
        Aminoacid/ModificationInfo.h
        Aminoacid/ModificationInfo.cpp
        ArgParse/Args.cpp
        ArgParse/Args.h
        Aminoacid/MonomerInfo.h
        Aminoacid/MonomerInfo.cpp
        NRP/MonomericNRPBuilder.h
        NRP/MonomericNRPBuilder.cpp)

add_library(nrps_matcher_lib ${SOURCE_FILES})
add_executable(NRPsMatcher main.cpp)
add_executable(ParseGraph parse_graphs.cpp)

target_link_libraries(NRPsMatcher logger_lib nrps_matcher_lib)
target_link_libraries(ParseGraph logger_lib nrps_matcher_lib)

install(TARGETS NRPsMatcher
        DESTINATION bin
        COMPONENT runtime)