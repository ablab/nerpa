//
// Created by olga on 26.03.19.
//

#include "ScoreFullMatch.h"

bool matcher::ScoreFullMatch::getScoreForSegment(const std::vector<aminoacid::Aminoacid> &amns,
                                                 const nrpsprediction::NRPsPrediction& prediction, int part_id,
                                                 double &score) const {
    nrpsprediction::NRPsPart part = prediction.getNrpsParts()[part_id];
    std::vector<nrpsprediction::AminoacidPrediction> aminoacid_predictions = part.getAminoacidsPrediction();
    double segscor = 0;
    for (int j = 0; j < (int)aminoacid_predictions.size(); ++j) {
        if (!aminoacid_predictions[j].contain(amns[j])) {
            return false;
        }
        double  cur_sc = aaScore(aminoacid_predictions[j], amns[j]);

        segscor += cur_sc;
    }
    score = segscor;
    return true;
}
