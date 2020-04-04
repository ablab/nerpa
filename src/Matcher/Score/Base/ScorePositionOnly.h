//
// Created by olga on 16.03.19.
//

#ifndef NRPSMATCHER_SCOREPOSITIONONLY_H
#define NRPSMATCHER_SCOREPOSITIONONLY_H

#include "Matcher/Score/Base/Score.h"

namespace matcher {
    class ScorePositionOnly : public Score {
    public:
        double
        aaScore(const nrpsprediction::AAdomain_Prediction &apred, const aminoacid::Aminoacid &aminoacid) const override;
    };
}

#endif //NRPSMATCHER_SCOREPOSITIONONLY_H
