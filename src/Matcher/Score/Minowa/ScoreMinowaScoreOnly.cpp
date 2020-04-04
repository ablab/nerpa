//
// Created by olga on 18.03.19.
//

#include "ScoreMinowaScoreOnly.h"

namespace matcher {
    double matcher::ScoreMinowaScoreOnly::aaScore(const nrpsprediction::AAdomain_Prediction &apred,
                                         const aminoacid::Aminoacid &aminoacid) const {
        nrpsprediction::AAdomain_Prediction::AminoacidProb prob = apred.getAminoacid(aminoacid);

        return std::min(1., prob.prob/350.);
    }
}