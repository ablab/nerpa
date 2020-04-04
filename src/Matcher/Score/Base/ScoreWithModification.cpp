//
// Created by olga on 01.02.19.
//

#include <cmath>
#include "ScoreWithModification.h"

namespace matcher {
    bool ScoreWithModification::getScoreForSegment(const std::vector<aminoacid::Aminoacid> &amns,
                                                   const nrpsprediction::BgcPrediction& prediction, int part_id,
                                                   double &score) const {
        nrpsprediction::OrfPrediction part = prediction.getOrfs()[part_id];
        std::vector<nrpsprediction::AAdomainPrediction> aminoacid_predictions = part.getAAdomainPrediction();
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

    double ScoreWithModification::aaScore(const nrpsprediction::AAdomainPrediction &apred,
                                          const aminoacid::Aminoacid &aminoacid) const {
        nrpsprediction::AAdomainPrediction::AminoacidProb probRes;
        std::pair<int, int> posRes;
        return getTheBestAAInPred(apred, aminoacid, probRes, posRes).first;
    }

    bool ScoreWithModification::getScore(const aminoacid::Aminoacid &nrpAA, const aminoacid::Aminoacid &predAA,
                                           const nrpsprediction::AAdomainPrediction::AminoacidProb &prob,
                                           const std::pair<int, int> &pos,
                                           double& score) const {
        aminoacid::Formula formula = (nrpAA - predAA);
        aminoacid::Modification modification(formula);
        if (modification.getId() == aminoacid::ModificationInfo::MODIFICATION_CNT) {
            return false;
        } else {
            double modCoeff = modification.getScore(predAA.get_id());
            if (modCoeff < 0) {
                return false;
            }
            if (baseScore->getScore(nrpAA, predAA, prob, pos, score)) {
                score *= modCoeff;
                return true;
            } else {
                return false;
            }
        }
    }

    std::pair<double, aminoacid::Aminoacid>
    ScoreWithModification::getTheBestAAInPred(const nrpsprediction::AAdomainPrediction &apred,
                                              const aminoacid::Aminoacid &aminoacid,
                                              nrpsprediction::AAdomainPrediction::AminoacidProb &probRes,
                                              std::pair<int, int> &posRes) const {
        auto AAprobs = apred.getAAPrediction();
        aminoacid::Aminoacid theBest;
        int bg = 0;
        int ed = 0;
        double maxScore = -1;
        bool has_match = false;
        for (int i = 0; i < AAprobs.size(); ++i) {
            if (i == 0 || fabs(AAprobs[i].prob - AAprobs[bg].prob) > EPS) {
                bg = i;
                ed = i;
                while (ed < AAprobs.size() && fabs(AAprobs[ed].prob - AAprobs[bg].prob) < EPS) {
                    ed += 1;
                }
            }

            double curScore;
            bool is_match = getScore(aminoacid, AAprobs[i].aminoacid, AAprobs[i], std::make_pair(bg, ed - 1), curScore);
            if (is_match && maxScore < curScore) {
                has_match = true;
                maxScore = curScore;
                theBest = AAprobs[i].aminoacid;
                probRes = AAprobs[i];
                aminoacid::Modification mod(aminoacid - AAprobs[i].aminoacid);
                if (aminoacid::ModificationInfo::NAMES[mod.getId()] != "empty") {
                    theBest.addModification(mod);
                }
                auto cur = std::make_pair(bg, ed - 1);
                posRes.swap(cur);
            }
        }

        if (has_match == false) {
            posRes = std::make_pair(-1, -1);
            return std::make_pair(Mismatch(), aminoacid::Aminoacid("none"));
        }

        return std::make_pair(maxScore, theBest);
    }

    ScoreWithModification::ScoreWithModification(std::unique_ptr<Score> base) : Score(std::move(base)) {}
}
