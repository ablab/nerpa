//
// Created by olga on 26.03.19.
//

#include "ScoreFullMatch.h"

bool matcher::ScoreFullMatch::getScoreForSegment(const std::vector<aminoacid::Aminoacid> &amns,
                                                 const nrpsprediction::BGC_Prediction& prediction, int part_id,
                                                 double &score) const {
    nrpsprediction::ORF_Prediction part = prediction.getNrpsParts()[part_id];
    std::vector<nrpsprediction::AAdomain_Prediction> aminoacid_predictions = part.getAminoacidsPrediction();
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
