//
// Created by olga on 22.01.19.
//

#ifndef NRPSMATCHER_SCORE_H
#define NRPSMATCHER_SCORE_H

#include <NRP/NRP.h>
#include <Matcher/Segment.h>

namespace matcher {
    class Score {
    public:
        double minScore(const int len) const {
            return -len - 1;
        }

        double openGap() const {
            return -1;
        }

        double continueGap() const {
            return 0;
        }

        double addSegment(Segment seg) const {
            return seg.scor - 1;
        }

        bool getScoreForSegment(const std::vector<aminoacid::Aminoacids::Aminoacid>& amns,
                                const nrpsprediction::NRPsPart& part, double& score) const;

        double aaScore(const nrpsprediction::AminoacidPrediction &apred,
                       const aminoacid::Aminoacids::Aminoacid &aminoacid) const;
    };
}


#endif //NRPSMATCHER_SCORE_H
