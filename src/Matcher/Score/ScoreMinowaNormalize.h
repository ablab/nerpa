#ifndef NRPSMATCHER_SCOREMINOWANORMALIZE_H
#define NRPSMATCHER_SCOREMINOWANORMALIZE_H

#include "Score.h"
#include "ScoreMinowa.h"

namespace matcher {
    class ScoreMinowaNormalize : public ScoreMinowa {
    public:
        double resultScore(double score, const int len,
                           const std::vector<Segment>& matched_parts,
                           const nrpsprediction::NRPsPrediction& prediction,
                           const nrp::NRP& nrp) const override;

    };
}


#endif //NRPSMATCHER_SCOREMINOWANORMALIZE_H
