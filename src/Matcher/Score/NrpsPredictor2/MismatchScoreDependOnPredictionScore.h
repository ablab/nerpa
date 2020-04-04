#ifndef NERPA_MISMATCHSCOREDEPENDONPREDICTIONSCORE_H
#define NERPA_MISMATCHSCOREDEPENDONPREDICTIONSCORE_H

#include <Matcher/Score/Base/Score.h>

namespace matcher {
    class MismatchScoreDependOnPredictionScore : public Score {
    public:
        explicit MismatchScoreDependOnPredictionScore(std::unique_ptr<Score> base) : Score(std::move(base)) {}

        double Mismatch(const aminoacid::Aminoacid &structure_aa,
                        const nrpsprediction::AAdomainPrediction &aa_prediction) const override;

    };

    double MismatchScoreDependOnPredictionScore::Mismatch(const aminoacid::Aminoacid &structure_aa,
                                           const nrpsprediction::AAdomainPrediction &aa_prediction) const {
        if (aa_prediction.getAAPrediction().size() == 0) {
            return 0;
        }

        double mismatch_score[11] = {0., 0., 0., 0., 0., 0., 0., -0.1, -0.25, -0.5, -1};
        return mismatch_score[int(aa_prediction.getAAPrediction()[0].prob/10)];
    }
}

#endif //NERPA_MISMATCHSCOREDEPENDONPREDICTIONSCORE_H
