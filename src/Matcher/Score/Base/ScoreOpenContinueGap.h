//
// Created by olga on 13.07.19.
//

#ifndef NRPSMATCHER_SCOREOPENCONTINUEGAP_H
#define NRPSMATCHER_SCOREOPENCONTINUEGAP_H

#include "Score.h"

namespace matcher {
    class ScoreOpenContinueGap : public Score {
    protected:
        double openGapScore = 0;
        double continueGapScore = 0;
    public:
        ScoreOpenContinueGap(double openGap, double continueGap, std::unique_ptr<Score> base) : Score(std::move(base)) {
            openGapScore = openGap;
            continueGapScore = continueGap;
        }

        ScoreOpenContinueGap(std::unique_ptr<Score> base) : Score(std::move(base)) {}

        double openGap() const override {
            return openGapScore;
        }

        double continueGap() const override {
            return continueGapScore;
        }
    };
}


#endif //NRPSMATCHER_SCOREOPENCONTINUEGAP_H
