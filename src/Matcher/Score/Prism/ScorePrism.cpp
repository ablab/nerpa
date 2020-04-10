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

    double ScorePrism::Mismatch(const aminoacid::Aminoacid &structure_aa,
                                const nrpsprediction::AAdomainPrediction &aa_prediction) const {
        if (baseScore != nullptr) {
            return baseScore->Mismatch(structure_aa, aa_prediction);
        } else {
            double mismatch_score[11];
            mismatch_score[10] = mismatch;
            for (int i = 9; i >= 0; --i) {
                mismatch_score[i] = mismatch_score[i + 1]/2;
            }

            if (int(aa_prediction.getAAPrediction()[0].prob/50) > 10) {
                return mismatch;
            }

            return mismatch_score[int(aa_prediction.getAAPrediction()[0].prob/50)];
        }
    }

    ScorePrism::ScorePrism(double mismatch) : Score(mismatch) {}
}