//
// Created by olga on 04.07.19.
//

#include "InDelMatcher.h"

namespace  matcher {
    matcher::MatcherBase::Match
    matcher::InDelMatcher::getMatch(const nrp::NRP *nrp, const nrpsprediction::NRPsPrediction *prediction,
                                    const matcher::Score *score) {
        MatcherBase::Match theBestMatch = innerMatcher->getMatch(nrp, prediction, score);

        return theBestMatch;
    }
}
