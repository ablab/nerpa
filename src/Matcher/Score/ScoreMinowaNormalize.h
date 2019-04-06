#ifndef NRPSMATCHER_SCOREMINOWANORMALIZE_H
#define NRPSMATCHER_SCOREMINOWANORMALIZE_H

#include "Score.h"
#include "ScoreMinowa.h"

namespace matcher {
    class ScoreMinowaNormalize : public ScoreMinowa {
    public:
        double resultScore(double score, const int len) const override;

    };
}


#endif //NRPSMATCHER_SCOREMINOWANORMALIZE_H
