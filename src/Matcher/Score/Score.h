//
// Created by olga on 22.01.19.
//

#ifndef NRPSMATCHER_SCORE_H
#define NRPSMATCHER_SCORE_H

#include <NRP/NRP.h>
#include <Matcher/Segment.h>

namespace matcher {
    class Score {
    public:
        double minScore(const nrp::NRP& nrp) const {
            return -nrp.getLen() - 1;
        }

        double openGap() const {
            return -1;
        }

        double continueGap() const {
            return 0;
        }

        double addSegment(Segment seg) const {
            return seg.scor - 1;
        }

        bool getScoreForSegment(const std::vector<aminoacid::Aminoacids::Aminoacid>& amns,
                                const nrpsprediction::NRPsPart& part, double& score) const {
            std::vector<nrpsprediction::AminoacidPrediction> aminoacid_predictions = part.getAminoacidsPrediction();
            int cnt_mismatch = 0;
            int g = 0;
            double segscor = 0;
            for (int j = 0; j < (int)aminoacid_predictions.size() && cnt_mismatch < 2; ++j) {
                if (!aminoacid_predictions[j].contain(amns[j])) {
                    cnt_mismatch += 1;
                }
                segscor += aminoacid_predictions[j].getScore(amns[j]);
            }

            if (cnt_mismatch == 0 || (cnt_mismatch == 1 && aminoacid_predictions.size() > 4)) {
                score = segscor;
                return true;
            }

            return false;
        }
    };
}


#endif //NRPSMATCHER_SCORE_H
