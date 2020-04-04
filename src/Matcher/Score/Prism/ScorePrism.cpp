//
// Created by olga on 02.03.19.
//

#include "ScorePrism.h"

namespace matcher {
    bool ScorePrism::getScore(const aminoacid::Aminoacid &nrpAA, const aminoacid::Aminoacid &predAA,
                                const nrpsprediction::AAdomainPrediction::AminoacidProb &prob,
                                const std::pair<int, int> &pos,
                                double& score) const {
        if (pos.first == -1) {
            return false;
        } else {
            score = prob.prob/500.;
        }
        return true;
    }
}