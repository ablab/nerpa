//
// Created by olga on 01.02.19.
//

#ifndef NRPSMATCHER_SCOREWITHMODIFICATION_H
#define NRPSMATCHER_SCOREWITHMODIFICATION_H

#include "Matcher/Score/Base/Score.h"

namespace matcher {
    class ScoreWithModification : public Score {
    public:
        ScoreWithModification(std::unique_ptr<Score> base);

    public:
        bool getScoreForSegment(const std::vector<aminoacid::Aminoacid> &amns,
                                const nrpsprediction::BGC_Prediction& prediction, int part_id,
                                double &score) const override;

        double
        aaScore(const nrpsprediction::AAdomain_Prediction &apred, const aminoacid::Aminoacid &aminoacid) const override;

        std::pair<double, aminoacid::Aminoacid> getTheBestAAInPred(const nrpsprediction::AAdomain_Prediction &apred,
                                                                   const aminoacid::Aminoacid &aminoacid,
                                                                   nrpsprediction::AAdomain_Prediction::AminoacidProb &probRes,
                                                                   std::pair<int, int> &posRes) const override;

        virtual double getScore(const aminoacid::Aminoacid& nrpAA,
                        const aminoacid::Aminoacid& predAA,
                        const nrpsprediction::AAdomain_Prediction::AminoacidProb& prob,
                        const std::pair<int, int>& pos) const;

        double EPS = 1e-5;
    };
}


#endif //NRPSMATCHER_SCOREWITHMODIFICATION_H
