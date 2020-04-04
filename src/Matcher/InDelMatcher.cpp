//
// Created by olga on 04.07.19.
//

#include <NRP/NRPLine.h>
#include <NRP/NRPCycle.h>
#include "InDelMatcher.h"

namespace  matcher {
    matcher::MatcherBase::Match
    matcher::InDelMatcher::getMatch(std::shared_ptr<nrp::NRP> nrp, const nrpsprediction::BGC_Prediction *prediction,
                                    const matcher::Score *score) {
        MatcherBase::Match theBestMatch = innerMatcher->getMatch(nrp, prediction, score);
        if (deletion) {
            MatcherBase::Match delBestMatch = getDeleteMatch(nrp, prediction, score);
            delBestMatch.setScore(score->InDelScore(delBestMatch.score(), nrp->getLen() - 1));
            if (delBestMatch.score() > theBestMatch.score()) {
                theBestMatch = delBestMatch;
            }
        }
        if (insertion) {
            MatcherBase::Match inBestMatch = getInsertMatch(nrp, prediction, score);
            inBestMatch.setScore(score->InDelScore(inBestMatch.score(), nrp->getLen() - 1));
            if (inBestMatch.score() > theBestMatch.score()) {
                theBestMatch = inBestMatch;
            }
        }

        return theBestMatch;
    }

    MatcherBase::Match
    InDelMatcher::getDeleteMatch(std::shared_ptr<nrp::NRP> nrp, const nrpsprediction::BGC_Prediction *prediction,
                                 const matcher::Score *score) {
        if (nrp->getType() == nrp::NRP::cycle) {
            MatcherBase::Match res = innerMatcher->getMatch(nrp, prediction, score);
            int len = nrp->getLen();
            for (int i = 0; i < len; ++i) {
                std::shared_ptr<nrp::NRP> nrpCycle = std::make_shared<nrp::NRPCycle>(*nrp);
                nrpCycle->deleteAA(i);
                MatcherBase::Match match = innerMatcher->getMatch(nrpCycle, prediction, score);
                if (match.score() > res.score()) {
                    res = match;
                }
            }
            return res;
        } else if (nrp->getType() == nrp::NRP::line) {
            MatcherBase::Match res = innerMatcher->getMatch(nrp, prediction, score);
            int len = nrp->getLen();
            for (int i = 0; i < len; ++i) {
                std::shared_ptr<nrp::NRP> nrpLine = std::make_shared<nrp::NRPLine>(*nrp);
                nrpLine->deleteAA(i);
                MatcherBase::Match match = innerMatcher->getMatch(nrpLine, prediction, score);
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

    MatcherBase::Match
    InDelMatcher::getInsertMatch(std::shared_ptr<nrp::NRP> nrp, const nrpsprediction::BGC_Prediction *prediction,
                                 const matcher::Score *score) {
        if (nrp->getType() == nrp::NRP::cycle) {
            MatcherBase::Match res = innerMatcher->getMatch(nrp, prediction, score);
            int len = nrp->getLen();
            for (int i = 0; i <= len; ++i) {
                std::shared_ptr<nrp::NRP> nrpCycle = std::make_shared<nrp::NRPCycle>(*nrp);
                nrpCycle->insertAA(i);
                MatcherBase::Match match = innerMatcher->getMatch(nrpCycle, prediction, score);
                if (match.score() > res.score()) {
                    res = match;
                }
            }
            return res;
        } else if (nrp->getType() == nrp::NRP::line) {
            MatcherBase::Match res = innerMatcher->getMatch(nrp, prediction, score);
            int len = nrp->getLen();
            for (int i = 0; i < len; ++i) {
                std::shared_ptr<nrp::NRP> nrpLine = std::make_shared<nrp::NRPLine>(*nrp);
                nrpLine->insertAA(i);
                MatcherBase::Match match = innerMatcher->getMatch(nrpLine, prediction, score);
                if (match.score() > res.score()) {
                    res = match;
                }
            }
            return res;
        } else {
            MatcherBase::Match firstMatch = getInsertMatch(nrp->getLines()[0], prediction, score);
            MatcherBase::Match secondMatch = getInsertMatch(nrp->getLines()[1], prediction, score);
            if (firstMatch.score() > secondMatch.score()) {
                return firstMatch;
            } else {
                return secondMatch;
            }
        }
    }
}
