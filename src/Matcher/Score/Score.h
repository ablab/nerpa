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
        Score();
        virtual double minScore(const int len) const {
            return -len - 1;
        }

        virtual double openGap() const {
            return -1;
        }

        virtual double continueGap() const {
            return 0;
        }

        virtual double addSegment(Segment seg) const {
            return seg.scor - 1;
        }

        virtual bool getScoreForSegment(const std::vector<aminoacid::Aminoacid>& amns,
                                const nrpsprediction::NRPsPart& part, double& score) const;

        virtual double aaScore(const nrpsprediction::AminoacidPrediction &apred,
                       const aminoacid::Aminoacid &aminoacid) const;

        virtual std::pair<double, aminoacid::Aminoacid> getTheBestAAInPred(const nrpsprediction::AminoacidPrediction &apred,
                                                                           const aminoacid::Aminoacid &aminoacid,
                                                                           nrpsprediction::AminoacidPrediction::AminoacidProb &probRes,
                                                                           std::pair<int, int> &posRes) const;

    protected:
        double posscore[100];
    };
}


#endif //NRPSMATCHER_SCORE_H
