//
// Created by olga on 28.02.19.
//

#include "ScoreMinowa.h"

double matcher::ScoreMinowa::aaScore(const nrpsprediction::AminoacidPrediction &apred,
                                     const aminoacid::Aminoacid &aminoacid) const {
    std::pair<int, int> position = apred.getAmnAcidPos(aminoacid);
    nrpsprediction::AminoacidPrediction::AminoacidProb prob = apred.getAminoacid(aminoacid);

    if (position.first == -1) {
        return -1;
    } else {
        int mdpos = (position.first + position.second)/2;
        if (mdpos >= 10) {
            return 0;
        }
        return prob.prob/400. * posscore[mdpos];
    }
}
