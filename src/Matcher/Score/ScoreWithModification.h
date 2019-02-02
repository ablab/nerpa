//
// Created by olga on 01.02.19.
//

#ifndef NRPSMATCHER_SCOREWITHMODIFICATION_H
#define NRPSMATCHER_SCOREWITHMODIFICATION_H

#include "Score.h"

namespace matcher {
    class ScoreWithModification : public Score {
    public:
        double minScore(const int len) const override;

        double openGap() const override;

        double continueGap() const override;

        double addSegment(Segment seg) const override;

        bool getScoreForSegment(const std::vector<aminoacid::Aminoacid> &amns, const nrpsprediction::NRPsPart &part,
                                double &score) const override;

        double
        aaScore(const nrpsprediction::AminoacidPrediction &apred, const aminoacid::Aminoacid &aminoacid) const override;

    private:


    };
}


#endif //NRPSMATCHER_SCOREWITHMODIFICATION_H
