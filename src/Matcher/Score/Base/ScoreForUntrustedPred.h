//
// Created by olga on 07.12.19.
//

#ifndef NERPA_SCOREFORUNTRUSTEDPRED_H
#define NERPA_SCOREFORUNTRUSTEDPRED_H

#include "Matcher/Score/Base/Score.h"

namespace matcher {
    class ScoreForUntrustedPred : public Score {
    private:
        double prob_threshold = 300;
    public:
        ScoreForUntrustedPred(std::unique_ptr<Score> base) : Score(std::move(base)) {}
        ScoreForUntrustedPred(std::unique_ptr<Score> base, const std::string& predictor) : Score(std::move(base)) {
            if (predictor == "MINOWA") {
                prob_threshold = 250.;
            } else if (predictor == "NRPSPREDICTOR2") {
                prob_threshold = 70.;
            } else if (predictor == "SANDPUMA") {
                prob_threshold = 70.;
            } else if (predictor == "PRISM") {
                prob_threshold = 300.;
            }
        }

        double
        aaScore(const nrpsprediction::AAdomain_Prediction &apred, const aminoacid::Aminoacid &aminoacid) const override {
            if (aminoacid.get_possible_name() == "none") {
                return 0;
            }

            auto AAprobs = apred.getAAPrediction();
            assert(baseScore != nullptr);
            auto aa_score = baseScore->aaScore(apred, aminoacid);

            if (AAprobs.size() == 0 || AAprobs[0].prob < prob_threshold) {
                return std::max(0., aa_score);
            }

            return aa_score;
        }

        std::pair<double, aminoacid::Aminoacid> getTheBestAAInPred(const nrpsprediction::AAdomain_Prediction &apred,
                                                                   const aminoacid::Aminoacid &aminoacid,
                                                                   nrpsprediction::AAdomain_Prediction::AminoacidProb &probRes,
                                                                   std::pair<int, int> &posRes) const override {
            if (aminoacid.get_possible_name() == "none") {
                return std::make_pair(0., aminoacid::Aminoacid("none"));
            }

            auto AAprobs = apred.getAAPrediction();
            assert(baseScore != nullptr);
            auto best_aa = baseScore->getTheBestAAInPred(apred, aminoacid, probRes, posRes);

            if (AAprobs.size() == 0 || AAprobs[0].prob < prob_threshold) {
                if (best_aa.first < 0) {
                    return std::make_pair(0., aminoacid::Aminoacid("none"));
                }
            }

            return best_aa;
        };

        double EPS = 1e-5;
    };
}


#endif //NERPA_SCOREFORUNTRUSTEDPRED_H
