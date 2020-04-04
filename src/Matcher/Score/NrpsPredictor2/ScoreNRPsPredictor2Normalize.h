//
// Created by olga on 09.04.19.
//

#ifndef NRPSMATCHER_SCORENRPSPREDICTOR2NORMALIZE_H
#define NRPSMATCHER_SCORENRPSPREDICTOR2NORMALIZE_H

#include "Matcher/Score/Base/Score.h"

namespace matcher {
    class ScoreNRPsPredictor2Normalize : public Score {
        double resultScore(double score, const int len,
                           const std::vector<Segment>& matched_parts,
                           const nrpsprediction::BGC_Prediction& prediction,
                           const nrp::NRP& nrp) const override {
            return score / maxScore(len);
        }

        double openGap() const override {
            return 0;
        }

        double continueGap() const override {
            return 0;//(double)1/5; //-(double)1/10;
        }
    };
}


#endif //NRPSMATCHER_SCORENRPSPREDICTOR2NORMALIZE_H
