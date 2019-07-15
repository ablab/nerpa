//
// Created by olga on 13.07.19.
//

#ifndef NRPSMATCHER_SCORESINGLEUNIT_H
#define NRPSMATCHER_SCORESINGLEUNIT_H

#include "Score.h"

namespace matcher {
    class ScoreSingleUnit : public Score {
    public:
        ScoreSingleUnit(std::unique_ptr<Score> base) : Score(std::move(base)) {}

    public:
        double singleUnitScore(const nrpsprediction::AminoacidPrediction &apred,
                               const aminoacid::Aminoacid &aminoacid) const override {
            return baseScore->aaScore(apred, aminoacid) * 0.5;
        }
    };
}


#endif //NRPSMATCHER_SCORESINGLEUNIT_H
