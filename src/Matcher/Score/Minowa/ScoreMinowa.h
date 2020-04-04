//
// Created by olga on 28.02.19.
//

#ifndef NRPSMATCHER_SCOREMINOWA_H
#define NRPSMATCHER_SCOREMINOWA_H

#include "Matcher/Score/Base/Score.h"

namespace matcher {
    class ScoreMinowa : public Score {
    public:
        double getScore(const aminoacid::Aminoacid &nrpAA, const aminoacid::Aminoacid &predAA,
                        const nrpsprediction::AAdomain_Prediction::AminoacidProb &prob,
                        const std::pair<int, int> &pos) const override;
    };
}


#endif //NRPSMATCHER_SCOREMINOWA_H
