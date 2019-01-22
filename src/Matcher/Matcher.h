//
// Created by olga on 22.01.19.
//

#ifndef NRPSMATCHER_MATCHER_H
#define NRPSMATCHER_MATCHER_H

#include <NRPsPrediction/AminoacidPrediction.h>
#include <NRP/NRP.h>

namespace matcher {
    class Matcher {
    private:
        const nrp::NRP& nrp;
        const nrpsprediction::NRPsPrediction& prediction;
    public:
        Matcher(const nrp::NRP &nrp, const nrpsprediction::NRPsPrediction& prediction):
                nrp(nrp), prediction(prediction) {}

        nrp::NRP::Match getMatch() const;
    private:
        nrp::NRP::Match getLineMatch() const;
        nrp::NRP::Match getCycleMatch() const;
        nrp::NRP::Match getBranchMatch() const;


        nrp::NRP::Match updateMatch(const nrpsprediction::NRPsPrediction& nrPsPrediction, nrp::NRP::Match match, int bg) const;
    };
}


#endif //NRPSMATCHER_MATCHER_H
