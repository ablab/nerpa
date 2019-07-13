//
// Created by olga on 28.02.19.
//

#ifndef NRPSMATCHER_SCOREMINOWA_H
#define NRPSMATCHER_SCOREMINOWA_H

#include "Matcher/Score/Base/Score.h"
#include "ScoreMinowaWithModification.h"

namespace matcher {
    class ScoreMinowa : public ScoreMinowaWithModification {
    public:
        double
        aaScore(const nrpsprediction::AminoacidPrediction &apred, const aminoacid::Aminoacid &aminoacid) const override;
    };
}


#endif //NRPSMATCHER_SCOREMINOWA_H
