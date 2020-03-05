//
// Created by olga on 22.01.19.
//

#include <Logger/logger.hpp>
#include "Score.h"

bool matcher::Score::getScoreForSegment(const std::vector<aminoacid::Aminoacid>& amns,
                                        const nrpsprediction::NRPsPrediction& prediction, int part_id,
                                        double& score) const {
    if (baseScore != nullptr) {
        return baseScore->getScoreForSegment(amns, prediction, part_id, score);
    } else {
        nrpsprediction::NRPsPart part = prediction.getNrpsParts()[part_id];
        std::vector<nrpsprediction::AminoacidPrediction> aminoacid_predictions = part.getAminoacidsPrediction();
        int cnt_mismatch = 0;
        int g = 0;
        double segscor = 0;
        for (int j = 0; j < (int) aminoacid_predictions.size() && cnt_mismatch < 2; ++j) {
            if (!aminoacid_predictions[j].contain(amns[j])) {
                cnt_mismatch += 1;
            }
            double cur_sc = aaScore(aminoacid_predictions[j], amns[j]);

            segscor += cur_sc;
        }

        if (cnt_mismatch == 0 || (cnt_mismatch == 1 && aminoacid_predictions.size() > 4)) {
            score = segscor;
            return true;
        }

        return false;
    }
}

double matcher::Score::aaScore(const nrpsprediction::AminoacidPrediction &apred,
                               const aminoacid::Aminoacid &aminoacid) const {
    if (baseScore != nullptr) {
        return baseScore->aaScore(apred, aminoacid);
    } else {
        std::pair<int, int> position = apred.getAmnAcidPos(aminoacid);
        nrpsprediction::AminoacidPrediction::AminoacidProb prob = apred.getAminoacid(aminoacid);
        return getScore(aminoacid, aminoacid, prob, position);
    }
}

matcher::Score::Score() {
    baseScore = nullptr;
    double curscore = 1;
    for (int i = 0; i < 100; ++i) {
        posscore[i] = curscore;
        curscore /= 1.25;
    }
}

matcher::Score::Score(std::unique_ptr<Score> base) : Score() {
    baseScore = std::move(base);
}

std::pair<double, aminoacid::Aminoacid>
matcher::Score::getTheBestAAInPred(const nrpsprediction::AminoacidPrediction &apred,
                                   const aminoacid::Aminoacid &aminoacid,
                                   nrpsprediction::AminoacidPrediction::AminoacidProb &probRes,
                                   std::pair<int, int> &posRes) const {
    if (baseScore != nullptr) {
        return baseScore->getTheBestAAInPred(apred, aminoacid, probRes, posRes);
    } else {
        double score = aaScore(apred, aminoacid);
        auto curpos = apred.getAmnAcidPos(aminoacid);
        std::swap(posRes, curpos);
        probRes = apred.getAminoacid(aminoacid);
        return std::pair<double, aminoacid::Aminoacid>(score, probRes.aminoacid);
    }
}

double matcher::Score::singleUnitScore(const nrpsprediction::AminoacidPrediction &apred,
                                       const aminoacid::Aminoacid &aminoacid) const {
    if (baseScore != nullptr) {
        return baseScore->singleUnitScore(apred, aminoacid);
    } else {
        return 0;
    }
}

double matcher::Score::getScore(const aminoacid::Aminoacid &nrpAA, const aminoacid::Aminoacid &predAA,
                                const nrpsprediction::AminoacidPrediction::AminoacidProb &prob,
                                const std::pair<int, int> &pos) const {
    if (baseScore != nullptr) {
        return baseScore->getScore(nrpAA, predAA, prob, pos);
    } else {
        if (pos.first == -1) {
            return Mismatch();
        } else {
            int mdpos = (pos.first + pos.second) / 2;
            if (mdpos >= 10) {
                return 0;
            }
            return prob.prob / 100. * posscore[mdpos];
        }
    }
}
