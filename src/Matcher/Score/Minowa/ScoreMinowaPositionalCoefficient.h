//
// Created by olga on 27.04.19.
//

#ifndef NRPSMATCHER_SCOREMINOWAPOSITIONALCOEFFICIENT_H
#define NRPSMATCHER_SCOREMINOWAPOSITIONALCOEFFICIENT_H

#include <algorithm>
#include "Matcher/Score/Base/Score.h"
#include "ScoreMinowa.h"

namespace matcher {
    class ScoreMinowaPositionalCoefficient : public ScoreMinowa {
    private:
        double positional_coefficient = 2;

        double getScoreForSubseg(const Segment& seg, const nrp::NRP& nrp,
                                 const nrpsprediction::BGC_Prediction &prediction, int part_id) const {
            double score = 0;
            std::vector<aminoacid::Aminoacid> amns;
            if (seg.rev) {
                for (int i = seg.r; i != seg.l; i = (i - 1 + nrp.getLen()) % nrp.getLen()) {
                    amns.push_back(nrp.getAminoacid(i));
                }
                amns.push_back(nrp.getAminoacid(seg.l));
            } else {
                for (int i = seg.l; i != seg.r; i = (i + 1) % nrp.getLen()) {
                    amns.push_back(nrp.getAminoacid(i));
                }
                amns.push_back(nrp.getAminoacid(seg.r));
            }

            Score::getScoreForSegment(amns, prediction, part_id, score);
            return score;
        }
    public:
        double minScore(const int len) const override {
            return 0;
        }

        double maxScore(const int len) const override {
            return len;
        }

        double openGap() const override {
            return 0;
        }

        double continueGap() const override {
            return 0;
        }

        double addSegment(Segment seg) const override {
            return seg.scor;
        }

        bool getScoreForSegment(const std::vector<aminoacid::Aminoacid> &amns,
                                const nrpsprediction::BGC_Prediction &prediction, int part_id,
                                double &score) const override {
            bool isSeg = Score::getScoreForSegment(amns, prediction, part_id, score);
            if (!isSeg) {
                return isSeg;
            }
            auto parts = prediction.getNrpsParts();
            int len = parts[part_id].getAminoacidsPrediction().size();
            auto cmp = [](nrpsprediction::ORF_Prediction a , nrpsprediction::ORF_Prediction b) -> bool {
                return a.getAminoacidsPrediction().size() > b.getAminoacidsPrediction().size();
            };
            std::sort(parts.begin(), parts.end(), cmp);
            auto eq = [](nrpsprediction::ORF_Prediction a , nrpsprediction::ORF_Prediction b) -> bool {
                return a.getAminoacidsPrediction().size() == b.getAminoacidsPrediction().size();
            };
            parts.resize(std::unique(parts.begin(), parts.end(), eq) - parts.begin());

            int cnt_longer = 0;
            for (int i = 0; i < parts.size(); ++i) {
                if (parts[i].getAminoacidsPrediction().size() > len) {
                    ++cnt_longer;
                }
            }

            for (int i = 0; i < cnt_longer; ++i) {
                score /= positional_coefficient;
            }
            return isSeg;
        }

        double resultScore(double score, const int len,
                           const std::vector<Segment>& matched_parts,
                           const nrpsprediction::BGC_Prediction& prediction,
                           const nrp::NRP& nrp) const override;
    };

    double ScoreMinowaPositionalCoefficient::resultScore(double score, const int len,
                                                         const std::vector<Segment>& matched_parts,
                                                         const nrpsprediction::BGC_Prediction& prediction,
                                                         const nrp::NRP& nrp) const {
        auto parts = prediction.getNrpsParts();
        std::vector<nrpsprediction::ORF_Prediction> not_matched_parts;
        int mi = 0, pi = 0;
        while (mi < matched_parts.size() && pi < parts.size()) {
            while (mi < matched_parts.size() && matched_parts[mi].part_id < pi) {
                ++mi;
            }
            while (pi < parts.size() && (mi == matched_parts[mi].part_id || pi < matched_parts[mi].part_id)) {
                not_matched_parts.push_back(parts[pi]);
                ++pi;
            }
            while (pi < parts.size() && mi < matched_parts.size() && pi == matched_parts[mi].part_id) {
                ++pi;
                ++mi;
            }
        }

        double resScore = 0;
        auto cmp = [](nrpsprediction::ORF_Prediction a , nrpsprediction::ORF_Prediction b) -> bool {
            return a.getAminoacidsPrediction().size() > b.getAminoacidsPrediction().size();
        };
        std::sort(not_matched_parts.begin(), not_matched_parts.end(), cmp);
        auto eq = [](nrpsprediction::ORF_Prediction a , nrpsprediction::ORF_Prediction b) -> bool {
            return a.getAminoacidsPrediction().size() == b.getAminoacidsPrediction().size();
        };
        not_matched_parts.resize(std::unique(not_matched_parts.begin(), not_matched_parts.end(), eq) -
                                 not_matched_parts.begin());
        for (int i = 0; i < matched_parts.size(); ++i) {
            int cntMore = 0;
            while (cntMore < not_matched_parts.size() &&
                    not_matched_parts[cntMore].getAminoacidsPrediction().size() >
                            (matched_parts[i].r - matched_parts[i].l + nrp.getLen()) % nrp.getLen() + 1) {
                ++cntMore;
            }

            double curSegScore = getScoreForSubseg(matched_parts[i], nrp, prediction, matched_parts[i].part_id);
            for (int j = 0; j < cntMore; ++j) {
                curSegScore /= positional_coefficient;
            }
            resScore += curSegScore;
        }
        assert(resScore >= score - 0.0001);

        return resScore/maxScore(len);
    }
}


#endif //NRPSMATCHER_SCOREMINOWAPOSITIONALCOEFFICIENT_H
