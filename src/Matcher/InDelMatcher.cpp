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
        MatcherBase::Match inBestMatch = getInsertMatch(nrp, prediction, score);

        if (theBestMatch.score() > delBestMatch.score() && theBestMatch.score() > inBestMatch.score()) {
            return theBestMatch;
        }

        if (delBestMatch.score() > inBestMatch.score()) {
            return delBestMatch;
        } else {
            return inBestMatch;
        }
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
                nrp::NRPLine nrpLine(*nrp);
                nrpLine.deleteAA(i);
                MatcherBase::Match match = innerMatcher->getMatch(&nrpLine, prediction, score);
                if (match.score() > res.score()) {
                    res = match;
                }
            }
            return res;
        } else {
            MatcherBase::Match firstMatch = getDeleteMatch(nrp->getLines()[0].get(), prediction, score);
            MatcherBase::Match secondMatch = getDeleteMatch(nrp->getLines()[1].get(), prediction, score);
            if (firstMatch.score() > secondMatch.score()) {
                return firstMatch;
            } else {
                return secondMatch;
            }
        }
    }

    MatcherBase::Match
    InDelMatcher::getInsertMatch(const nrp::NRP *nrp, const nrpsprediction::NRPsPrediction *prediction,
                                 const matcher::Score *score) {
        if (nrp->getType() == nrp::NRP::cycle) {
            MatcherBase::Match res = innerMatcher->getMatch(nrp, prediction, score);
            int len = nrp->getLen();
            for (int i = 0; i <= len; ++i) {
                nrp::NRPCycle nrpCycle(*nrp);
                nrpCycle.insertAA(i);
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
                nrp::NRPLine nrpLine(*nrp);
                nrpLine.insertAA(i);
                MatcherBase::Match match = innerMatcher->getMatch(&nrpLine, prediction, score);
                if (match.score() > res.score()) {
                    res = match;
                }
            }
            return res;
        } else {
            MatcherBase::Match firstMatch = getInsertMatch(nrp->getLines()[0].get(), prediction, score);
            MatcherBase::Match secondMatch = getInsertMatch(nrp->getLines()[1].get(), prediction, score);
            if (firstMatch.score() > secondMatch.score()) {
                return firstMatch;
            } else {
                return secondMatch;
            }
        }
    }
}
