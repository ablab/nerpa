//
// Created by olga on 28.02.19.
//

#include "ScoreMinowa.h"

bool matcher::ScoreMinowa::getScore(const aminoacid::Aminoacid &nrpAA, const aminoacid::Aminoacid &predAA,
                                      const nrpsprediction::AAdomainPrediction::AminoacidProb &prob,
                                      const std::pair<int, int> &pos,
                                      double& score) const {
    if (pos.first == -1) {
        return false;
    } else {
        int mdpos = (pos.first + pos.second)/2;
        if (mdpos >= 10) {
            score = 0;
            return true;
        }
        score = std::min(1., prob.prob/350.) * posscore[mdpos];
    }
    return true;
}
