//
// Created by olga on 18.03.19.
//

#include "ScoreMinowaScoreOnly.h"

namespace matcher {
    double matcher::ScoreMinowaScoreOnly::aaScore(const nrpsprediction::AminoacidPrediction &apred,
                                         const aminoacid::Aminoacid &aminoacid) const {
        nrpsprediction::AminoacidPrediction::AminoacidProb prob = apred.getAminoacid(aminoacid);

        return std::min(1., prob.prob/350.);
    }
}