#ifndef NRPSMATCHER_SCOREFULLMATCH_H
#define NRPSMATCHER_SCOREFULLMATCH_H

#include "Score.h"

namespace matcher {
    class ScoreFullMatch : public Score {
    public:
        bool getScoreForSegment(const std::vector<aminoacid::Aminoacid> &amns, const nrpsprediction::NRPsPart &part,
                                double &score) const override;

    };
}


#endif //NRPSMATCHER_SCOREFULLMATCH_H
