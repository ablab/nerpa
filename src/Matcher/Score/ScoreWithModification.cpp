//
// Created by olga on 01.02.19.
//

#include <cmath>
#include "ScoreWithModification.h"

namespace matcher {

    double ScoreWithModification::minScore(const int len) const {
        return Score::minScore(len);
    }

    double ScoreWithModification::openGap() const {
        return Score::openGap();
    }

    double ScoreWithModification::continueGap() const {
        return Score::continueGap();
    }

    double ScoreWithModification::addSegment(Segment seg) const {
        return Score::addSegment(seg);
    }

    bool ScoreWithModification::getScoreForSegment(const std::vector<aminoacid::Aminoacid> &amns,
                                                   const nrpsprediction::NRPsPrediction& prediction, int part_id,
                                                   double &score) const {
        nrpsprediction::NRPsPart part = prediction.getNrpsParts()[part_id];
        std::vector<nrpsprediction::AminoacidPrediction> aminoacid_predictions = part.getAminoacidsPrediction();
        int cnt_mismatch = 0;
        int g = 0;
        double segscor = 0;
        for (int j = 0; j < (int)aminoacid_predictions.size() && cnt_mismatch < 2; ++j) {
            double cur_sc = aaScore(aminoacid_predictions[j], amns[j]);
            if (cur_sc < 0) {
                cnt_mismatch += 1;
            }

            segscor += cur_sc;
        }
        if (cnt_mismatch == 0 || (cnt_mismatch == 1 && aminoacid_predictions.size() > 4)) {
            score = segscor;
            return true;
        }

        return false;
    }

    double ScoreWithModification::aaScore(const nrpsprediction::AminoacidPrediction &apred,
                                          const aminoacid::Aminoacid &aminoacid) const {
        nrpsprediction::AminoacidPrediction::AminoacidProb probRes;
        std::pair<int, int> posRes;
        return getTheBestAAInPred(apred, aminoacid, probRes, posRes).first;
    }

    double ScoreWithModification::getScore(const aminoacid::Aminoacid &nrpAA, const aminoacid::Aminoacid &predAA,
                                           const nrpsprediction::AminoacidPrediction::AminoacidProb &prob,
                                           const std::pair<int, int> &pos) const {
        aminoacid::Formula formula = (nrpAA - predAA);
        aminoacid::Modification modification(formula);
        if (modification.getId() == aminoacid::Modification::MODIFICATION_CNT) {
            return -1;
        }
        if (modification.getId() == aminoacid::Modification::empty ||
                (predAA.get_name() == "asn" && modification.getId() == aminoacid::Modification::OH) ||
                (predAA.get_name() == "glu" && modification.getId() == aminoacid::Modification::me3)) {
            int mdpos = (pos.first + pos.second)/2;
            if (mdpos >= 10) {
                return 0;
            }
            return prob.prob/100. * posscore[mdpos];
        }
    }

    std::pair<double, aminoacid::Aminoacid>
    ScoreWithModification::getTheBestAAInPred(const nrpsprediction::AminoacidPrediction &apred,
                                              const aminoacid::Aminoacid &aminoacid,
                                              nrpsprediction::AminoacidPrediction::AminoacidProb &probRes,
                                              std::pair<int, int> &posRes) const {
        auto AAprobs = apred.getAAPrediction();
        aminoacid::Aminoacid theBest;
        int bg = 0;
        int ed = 0;
        double maxScore = -1;
        for (int i = 0; i < AAprobs.size(); ++i) {
            if (i == 0 || fabs(AAprobs[i].prob - AAprobs[bg].prob) > EPS) {
                bg = i;
                ed = i;
                while (ed < AAprobs.size() && fabs(AAprobs[ed].prob - AAprobs[bg].prob) < EPS) {
                    ed += 1;
                }
            }
            double curScore = getScore(aminoacid, AAprobs[i].aminoacid, AAprobs[i], std::make_pair(bg, ed - 1));
            if (maxScore < curScore) {
                maxScore = curScore;
                theBest = AAprobs[i].aminoacid;
                probRes = AAprobs[i];
                aminoacid::Modification mod(aminoacid - AAprobs[i].aminoacid);
                if (mod.getId() != mod.empty) {
                    theBest.addModification(mod);
                }
                auto cur = std::make_pair(bg, ed - 1);
                posRes.swap(cur);
            }
        }

        return std::make_pair(maxScore, theBest);
    }
}
