//
// Created by olga on 02.03.19.
//

#include "ScorePrism.h"

namespace matcher {
    double
    ScorePrism::aaScore(const nrpsprediction::AminoacidPrediction &apred, const aminoacid::Aminoacid &aminoacid) const {
        std::pair<int, int> position = apred.getAmnAcidPos(aminoacid);
        nrpsprediction::AminoacidPrediction::AminoacidProb prob = apred.getAminoacid(aminoacid);

        if (position.first == -1) {
            return -1;
        } else {
            return prob.prob/500.;
        }
    }
}