//
// Created by olga on 01.02.19.
//

#include "ScoreWithModification.h"

namespace matcher {

    double ScoreWithModification::minScore(const int len) const {
        return Score::minScore(len);
    }

    double ScoreWithModification::openGap() const {
        return Score::openGap();
    }

    double ScoreWithModification::continueGap() const {
        return Score::continueGap();
    }

    double ScoreWithModification::addSegment(Segment seg) const {
        return Score::addSegment(seg);
    }

    bool ScoreWithModification::getScoreForSegment(const std::vector<aminoacid::Aminoacid> &amns,
                                                   const nrpsprediction::NRPsPart &part, double &score) const {
        return Score::getScoreForSegment(amns, part, score);
    }

    double ScoreWithModification::aaScore(const nrpsprediction::AminoacidPrediction &apred,
                                          const aminoacid::Aminoacid &aminoacid) const {
        return Score::aaScore(apred, aminoacid);
    }
}
