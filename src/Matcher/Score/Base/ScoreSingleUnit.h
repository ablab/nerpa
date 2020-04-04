//
// Created by olga on 13.07.19.
//

#ifndef NRPSMATCHER_SCORESINGLEUNIT_H
#define NRPSMATCHER_SCORESINGLEUNIT_H

#include "Score.h"

namespace matcher {
    class ScoreSingleUnit : public Score {
    private:
        double coeff = 0.1;
    public:
        ScoreSingleUnit(double coeff, std::unique_ptr<Score> base) : Score(std::move(base)) {
            this->coeff = coeff;
        }

        ScoreSingleUnit(std::unique_ptr<Score> base) : Score(std::move(base)) {}

    public:
        double singleUnitScore(const nrpsprediction::AAdomain_Prediction &apred,
                               const aminoacid::Aminoacid &aminoacid) const override {
            return baseScore->aaScore(apred, aminoacid) * coeff;
        }
    };
}


#endif //NRPSMATCHER_SCORESINGLEUNIT_H
