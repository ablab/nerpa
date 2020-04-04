//
// Created by olga on 13.07.19.
//

#ifndef NRPSMATCHER_SCORENORMALIZE_H
#define NRPSMATCHER_SCORENORMALIZE_H

#include "Score.h"

namespace matcher {
    class ScoreNormalize : public Score {
    public:
        ScoreNormalize(std::unique_ptr<Score> base);

        double resultScore(double score, const int len,
                           const std::vector<Segment>& matched_parts,
                           const nrpsprediction::BGC_Prediction& prediction,
                           const nrp::NRP& nrp) const override;
    };
}


#endif //NRPSMATCHER_SCORENORMALIZE_H
