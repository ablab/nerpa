//
// Created by olga on 18.03.19.
//

#ifndef NRPSMATCHER_SCOREMINOWASCOREONLY_H
#define NRPSMATCHER_SCOREMINOWASCOREONLY_H

#include "Score.h"

namespace matcher {
    class ScoreMinowaScoreOnly : public Score {
    public:
        double
        aaScore(const nrpsprediction::AminoacidPrediction &apred, const aminoacid::Aminoacid &aminoacid) const override;
    };
}

#endif //NRPSMATCHER_SCOREMINOWASCOREONLY_H
