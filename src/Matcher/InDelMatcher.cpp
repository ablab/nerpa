//
// Created by olga on 04.07.19.
//

#include <NRP/NRPLine.h>
#include <NRP/NRPCycle.h>
#include "InDelMatcher.h"

namespace  matcher {
    matcher::MatcherBase::Match
    matcher::InDelMatcher::getMatch(const nrp::NRP *nrp, const nrpsprediction::NRPsPrediction *prediction,
                                    const matcher::Score *score) {
        MatcherBase::Match theBestMatch = innerMatcher->getMatch(nrp, prediction, score);
        MatcherBase::Match delBestMatch = getDeleteMatch(nrp, prediction, score);

        if (delBestMatch.score() > theBestMatch.score()) {
            return delBestMatch;
        }

        return theBestMatch;
    }

    MatcherBase::Match
    InDelMatcher::getDeleteMatch(const nrp::NRP *nrp, const nrpsprediction::NRPsPrediction *prediction,
                                 const matcher::Score *score) {
        if (nrp->getType() == nrp::NRP::cycle) {
            MatcherBase::Match res = innerMatcher->getMatch(nrp, prediction, score);
            int len = nrp->getLen();
            for (int i = 0; i < len; ++i) {
                nrp::NRPCycle nrpCycle(*nrp);
                nrpCycle.deleteAA(i);
                MatcherBase::Match match = innerMatcher->getMatch(&nrpCycle, prediction, score);
                if (match.score() > res.score()) {
                    res = match;
                }
            }
            return res;
        } else if (nrp->getType() == nrp::NRP::line) {
            MatcherBase::Match res = innerMatcher->getMatch(nrp, prediction, score);
            int len = nrp->getLen();
            for (int i = 0; i < len; ++i) {
                nrp::NRPLine nrpCycle(*nrp);
                nrpCycle.deleteAA(i);
                MatcherBase::Match match = innerMatcher->getMatch(&nrpCycle, prediction, score);
                if (match.score() > res.score()) {
                    res = match;
                }
            }
            return res;
        } else {
            MatcherBase::Match firstMatch = getDeleteMatch(nrp->getLines()[0], prediction, score);
            MatcherBase::Match secondMatch = getDeleteMatch(nrp->getLines()[1], prediction, score);
            if (firstMatch.score() > secondMatch.score()) {
                return firstMatch;
            } else {
                return secondMatch;
            }
        }
    }
}
