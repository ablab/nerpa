//
// Created by olga on 28.02.19.
//

#ifndef NRPSMATCHER_SCOREMINOWA_H
#define NRPSMATCHER_SCOREMINOWA_H

#include "Score.h"

namespace matcher {
    class ScoreMinowa : public Score {
    public:
        double
        aaScore(const nrpsprediction::AminoacidPrediction &apred, const aminoacid::Aminoacid &aminoacid) const override;
    };
}


#endif //NRPSMATCHER_SCOREMINOWA_H
