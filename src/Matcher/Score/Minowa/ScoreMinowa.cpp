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

matcher::ScoreMinowa::ScoreMinowa(double mismatch) : Score(mismatch) {}

double matcher::ScoreMinowa::Mismatch(const aminoacid::Aminoacid &structure_aa,
                                      const nrpsprediction::AAdomainPrediction &aa_prediction) const {
    if (baseScore != nullptr) {
        return baseScore->Mismatch(structure_aa, aa_prediction);
    } else {
        double mismatch_score[11];
        mismatch_score[10] = mismatch;
        for (int i = 9; i >= 0; --i) {
            mismatch_score[i] = mismatch_score[i + 1]/2;
        }

        if (int(aa_prediction.getAAPrediction()[0].prob/35) > 10) {
            return mismatch;
        }

        return mismatch_score[int(aa_prediction.getAAPrediction()[0].prob/35)];
    }
}
