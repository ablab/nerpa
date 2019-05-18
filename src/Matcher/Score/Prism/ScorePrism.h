//
// Created by olga on 02.03.19.
//

#ifndef NRPSMATCHER_SCOREPRISM_H
#define NRPSMATCHER_SCOREPRISM_H

#include "Matcher/Score/Base/Score.h"

namespace matcher {
    class ScorePrism : public Score {
    public:
        double
        aaScore(const nrpsprediction::AminoacidPrediction &apred, const aminoacid::Aminoacid &aminoacid) const override;

    };
}


#endif //NRPSMATCHER_SCOREPRISM_H
