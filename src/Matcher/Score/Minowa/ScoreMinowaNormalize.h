#ifndef NRPSMATCHER_SCOREMINOWANORMALIZE_H
#define NRPSMATCHER_SCOREMINOWANORMALIZE_H

#include "Matcher/Score/Base/Score.h"
#include "ScoreMinowa.h"
#include "ScoreMinowaWithModification.h"

namespace matcher {
    class ScoreMinowaNormalize : public ScoreMinowaWithModification {
    public:
        double resultScore(double score, const int len,
                           const std::vector<Segment>& matched_parts,
                           const nrpsprediction::NRPsPrediction& prediction,
                           const nrp::NRP& nrp) const override;

    };
}


#endif //NRPSMATCHER_SCOREMINOWANORMALIZE_H
