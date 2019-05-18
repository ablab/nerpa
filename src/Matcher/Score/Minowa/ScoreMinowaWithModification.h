//
// Created by olga on 18.05.19.
//

#ifndef NRPSMATCHER_SCOREMINOWAWITHMODIFICATION_H
#define NRPSMATCHER_SCOREMINOWAWITHMODIFICATION_H


#include "Matcher/Score/Base/ScoreWithModification.h"

namespace matcher {
    class ScoreMinowaWithModification : public ScoreWithModification {
    protected:
        double getScore(const aminoacid::Aminoacid &nrpAA, const aminoacid::Aminoacid &predAA,
                        const nrpsprediction::AminoacidPrediction::AminoacidProb &prob,
                        const std::pair<int, int> &pos) const override;
    };
}


#endif //NRPSMATCHER_SCOREMINOWAWITHMODIFICATION_H
