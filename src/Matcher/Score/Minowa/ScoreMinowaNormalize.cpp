//
// Created by olga on 06.04.19.
//

#include "ScoreMinowaNormalize.h"

namespace matcher {
    double ScoreMinowaNormalize::resultScore(double score, const int len,
                                             const std::vector<Segment>& matched_parts,
                                             const nrpsprediction::NRPsPrediction& prediction,
                                             const nrp::NRP& nrp) const {
        return score / maxScore(len);
    }
};
