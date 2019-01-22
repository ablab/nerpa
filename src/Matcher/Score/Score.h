//
// Created by olga on 22.01.19.
//

#ifndef NRPSMATCHER_SCORE_H
#define NRPSMATCHER_SCORE_H

#include <NRP/NRP.h>

namespace matcher {
    class Score {
    public:
        double minScore(const nrp::NRP& nrp) const {
            return -nrp.getLen() - 1;
        }

        double openGap() const {
            return -1;
        }

        double continueGap() const {
            return 0;
        }

        double addSegment(nrp::NRP::Segment seg) const {
            return seg.scor - 1;
        }
    };
}


#endif //NRPSMATCHER_SCORE_H
