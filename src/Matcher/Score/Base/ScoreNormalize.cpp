//
// Created by olga on 13.07.19.
//

#include "ScoreNormalize.h"

namespace matcher {
    double ScoreNormalize::resultScore(double score, const int len,
                                       const std::vector<Segment>& matched_parts,
                                       const nrpsprediction::BgcPrediction& prediction,
                                       const nrp::NRP& nrp) const {
        return score / maxScore(len);
    }

    ScoreNormalize::ScoreNormalize(std::unique_ptr<Score> base) : Score(std::move(base)) {}
};
