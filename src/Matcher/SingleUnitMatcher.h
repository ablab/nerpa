//
// Created by olga on 21.07.19.
//

#ifndef NRPSMATCHER_SINGLEUNITMATCHER_H
#define NRPSMATCHER_SINGLEUNITMATCHER_H

#include <Aminoacid/Aminoacid.h>
#include "Matcher.h"

namespace matcher {
    typedef aminoacid::Aminoacid aacid;

    class SingleUnitMatcher : public Matcher {
    protected:
        void matchSingleUnits(Match& match, std::vector<bool>& used_pos) const;

        Match
        updateMatch(const nrpsprediction::NRPsPrediction &nrPsPrediction, matcher::MatcherBase::Match match, int bg,
                    std::vector<Segment> &matched_parts_id) const override;

    public:
        SingleUnitMatcher(const std::shared_ptr<nrp::NRP> &nrp, const nrpsprediction::NRPsPrediction *prediction,
                          const Score *score);

        SingleUnitMatcher();
    };
}


#endif //NRPSMATCHER_SINGLEUNITMATCHER_H
