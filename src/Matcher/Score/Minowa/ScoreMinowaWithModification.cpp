//
// Created by olga on 18.05.19.
//

#include "ScoreMinowaWithModification.h"

double matcher::ScoreMinowaWithModification::getScore(const aminoacid::Aminoacid &nrpAA, const aminoacid::Aminoacid &predAA,
                                      const nrpsprediction::AminoacidPrediction::AminoacidProb &prob,
                                      const std::pair<int, int> &pos) const {
    aminoacid::Formula formula = (nrpAA - predAA);
    aminoacid::Modification modification(formula);
    if (modification.getId() == aminoacid::Modification::MODIFICATION_CNT) {
        return -1;
    } else {
        int mdpos = (pos.first + pos.second)/2;
        if (mdpos >= 10) {
            return 0;
        }
        double modCoeff = 1;

        if (modification.getId() != aminoacid::Modification::empty && modification.getId() != aminoacid::Modification::methylation) {
            modCoeff = 0;
        }

        return std::min(1., prob.prob/350.) * posscore[mdpos] * modCoeff;
    }
}
