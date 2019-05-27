//
// Created by olga on 06.04.19.
//

#ifndef NRPSMATCHER_SCOREMINOWANORMALIZEWITHOUTAFFINGAP_H
#define NRPSMATCHER_SCOREMINOWANORMALIZEWITHOUTAFFINGAP_H

#include "Matcher/Score/Base/Score.h"
#include "ScoreMinowa.h"
#include "ScoreMinowaNormalize.h"

namespace matcher {
    class ScoreMinowaNormalizeWithoutAffinGap : public ScoreMinowaNormalize {
    public:
        double openGap() const override {
            return 0;
        }

        double continueGap() const override {
            return 0;//(double)1/5; //-(double)1/10;
        }

        /*double resultScore(double score, const int len) const {
            return score;
        }*/

        double singleUnitScore(const nrpsprediction::AminoacidPrediction &apred,
                               const aminoacid::Aminoacid &aminoacid) const override {
            return aaScore(apred, aminoacid) * 0.5;
        }
    };
}


#endif //NRPSMATCHER_SCOREMINOWANORMALIZEWITHOUTAFFINGAP_H
