//
// Created by olga on 02.03.19.
//

#include "ScorePrism.h"

namespace matcher {
    double ScorePrism::getScore(const aminoacid::Aminoacid &nrpAA, const aminoacid::Aminoacid &predAA,
                                const nrpsprediction::AminoacidPrediction::AminoacidProb &prob,
                                const std::pair<int, int> &pos) const {
        if (pos.first == -1) {
            return -1;
        } else {
            return prob.prob/500.;
        }
    }
}