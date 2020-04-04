//
// Created by olga on 28.02.19.
//

#ifndef NRPSMATCHER_SCOREMINOWA_H
#define NRPSMATCHER_SCOREMINOWA_H

#include "Matcher/Score/Base/Score.h"

namespace matcher {
    class ScoreMinowa : public Score {
    public:
        explicit ScoreMinowa(double mismatch);

        bool getScore(const aminoacid::Aminoacid &nrpAA, const aminoacid::Aminoacid &predAA,
                        const nrpsprediction::AAdomainPrediction::AminoacidProb &prob,
                        const std::pair<int, int> &pos,
                        double& score) const override;

        double Mismatch(const aminoacid::Aminoacid &structure_aa,
                        const nrpsprediction::AAdomainPrediction &aa_prediction) const override;
    };
}


#endif //NRPSMATCHER_SCOREMINOWA_H
