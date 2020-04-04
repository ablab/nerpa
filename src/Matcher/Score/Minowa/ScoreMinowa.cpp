//
// Created by olga on 28.02.19.
//

#include "ScoreMinowa.h"

double matcher::ScoreMinowa::getScore(const aminoacid::Aminoacid &nrpAA, const aminoacid::Aminoacid &predAA,
                                      const nrpsprediction::AAdomain_Prediction::AminoacidProb &prob,
                                      const std::pair<int, int> &pos) const {
    if (pos.first == -1) {
        return -1;
    } else {
        int mdpos = (pos.first + pos.second)/2;
        if (mdpos >= 10) {
            return 0;
        }
        return std::min(1., prob.prob/350.) * posscore[mdpos];
    }
}
