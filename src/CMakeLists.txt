include_directories(.)

add_subdirectory(Logger)
set(SOURCE_FILES main.cpp
        Aminoacid/Aminoacid.cpp
        Aminoacid/Aminoacid.h
        NRPsPrediction/AAdomain_Prediction.cpp
        NRPsPrediction/AAdomain_Prediction.h
        NRPsPrediction/ORF_Prediction.cpp
        NRPsPrediction/ORF_Prediction.h
        NRPsPrediction/BGC_Prediction.cpp
        NRPsPrediction/BGC_Prediction.h
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
        Matcher/OrderedGenesMatcher.h
        Matcher/OrderedGenesMatcher.cpp
        Matcher/MatcherBase.h
        Matcher/InDelMatcher.h
        Matcher/InDelMatcher.cpp
        Matcher/SingleUnitMatcher.cpp
        Matcher/SingleUnitMatcher.h
        Matcher/Score/Base/Score.cpp
        Matcher/Score/Base/Score.h
        Matcher/Segment.h
        NRPsPrediction/Builders/PredictionBuilderBase.h
        NRPsPrediction/Builders/PredictionBuilderBase.cpp
        NRPsPrediction/Builders/Nrpspredictor2Builder.h
        NRPsPrediction/Builders/Nrpspredictor2Builder.cpp
        NRPsPrediction/Builders/MinowaPredictionBuilder.h
        NRPsPrediction/Builders/MinowaPredictionBuilder.cpp
        NRPsPrediction/Builders/PrismPredictionBuilder.h
        NRPsPrediction/Builders/PrismPredictionBuilder.cpp
        NRPsPrediction/Builders/SandpumaPredictionBuilder.h
        NRPsPrediction/Builders/SandpumaPredictionBuilder.cpp
        Aminoacid/Formula.h
        Aminoacid/Formula.cpp
        Aminoacid/Modification.h
        Aminoacid/Modification.cpp
        Aminoacid/AminoacidInfo.h
        Aminoacid/AminoacidInfo.cpp
        Aminoacid/ModificationInfo.h
        Aminoacid/ModificationInfo.cpp
        Matcher/Score/Base/ScoreWithModification.cpp
        Matcher/Score/Base/ScoreWithModification.h
        Matcher/Score/Minowa/ScoreMinowa.cpp
        Matcher/Score/Minowa/ScoreMinowa.h
        Matcher/Score/Prism/ScorePrism.cpp
        Matcher/Score/Prism/ScorePrism.h
        Matcher/Score/Sandpuma/ScoreSandpuma.cpp
        Matcher/Score/Sandpuma/ScoreSandpuma.h
        Matcher/Score/Base/ScorePositionOnly.cpp
        Matcher/Score/Base/ScorePositionOnly.h
        Matcher/Score/Minowa/ScoreMinowaScoreOnly.cpp
        Matcher/Score/Minowa/ScoreMinowaScoreOnly.h
        Matcher/Score/Base/ScoreFullMatch.cpp
        Matcher/Score/Base/ScoreFullMatch.h
        Matcher/Score/NrpsPredictor2/ScoreNRPsPredictor2Normalize.h
        Matcher/Score/Minowa/ScoreMinowaPositionalCoefficient.h
        Matcher/Score/Base/ScoreNormalize.cpp
        Matcher/Score/Base/ScoreNormalize.h
        Matcher/Score/Base/ScoreOpenContinueGap.h
        Matcher/Score/Base/ScoreSingleUnit.h
        Matcher/Score/OrderedGenes/OrderedGenesScoreBase.h
        Matcher/Score/Base/ScoreForUntrustedPred.h
        Json/json.hpp
        ArgParse/Args.cpp
        ArgParse/Args.h)

add_library(nrps_matcher_lib ${SOURCE_FILES})
add_executable(NRPsMatcher main.cpp)
add_executable(ParseGraph parse_graphs.cpp)

target_link_libraries(NRPsMatcher logger_lib nrps_matcher_lib)
target_link_libraries(ParseGraph logger_lib nrps_matcher_lib)

install(TARGETS NRPsMatcher
        DESTINATION bin
        COMPONENT runtime)