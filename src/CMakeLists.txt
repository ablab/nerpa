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
        Matcher/MatcherBase.h
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
        Aminoacid/Modification.h
        Aminoacid/Modification.cpp
        Aminoacid/AminoacidInfo.h
        Aminoacid/AminoacidInfo.cpp
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
        Matcher/Score/Minowa/ScoreMinowaNormalize.h
        Matcher/Score/Minowa/ScoreMinowaNormalize.cpp
        Matcher/Score/Minowa/ScoreMinowaNormalizeWithoutAffinGap.h
        Matcher/Score/NrpsPredictor2/ScoreNRPsPredictor2Normalize.h
        Matcher/Score/Prism/ScorePrismNormalize.h
        Matcher/Score/Sandpuma/ScoreSandpumaNormalize.h
        Matcher/Score/Minowa/ScoreMinowaPositionalCoefficient.h
        Matcher/Score/Minowa/ScoreMinowaWithModification.h
        Matcher/Score/Minowa/ScoreMinowaWithModification.cpp
        Json/json.hpp)

add_library(nrps_matcher_lib ${SOURCE_FILES})
add_executable(NRPsMatcher main.cpp)

target_link_libraries(NRPsMatcher logger_lib nrps_matcher_lib)

install(TARGETS NRPsMatcher
        DESTINATION bin
        COMPONENT runtime)