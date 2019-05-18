//
// Created by olga on 10.04.19.
//

#ifndef NRPSMATCHER_SCORESANDPUMANORMALIZE_H
#define NRPSMATCHER_SCORESANDPUMANORMALIZE_H

#include "ScoreSandpuma.h"

namespace matcher {
    class ScoreSandpumaNormalize : public ScoreSandpuma {
        double resultScore(double score, const int len,
                           const std::vector<Segment>& matched_parts,
                           const nrpsprediction::NRPsPrediction& prediction,
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


#endif //NRPSMATCHER_SCORESANDPUMANORMALIZE_H
