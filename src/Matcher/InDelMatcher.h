//
// Created by olga on 04.07.19.
//

#ifndef NRPSMATCHER_INDELMATCHER_H
#define NRPSMATCHER_INDELMATCHER_H

#include "MatcherBase.h"
#include "Matcher.h"

namespace matcher {
    class InDelMatcher : public MatcherBase {
    private:
        MatcherBase* innerMatcher = nullptr;

        ~InDelMatcher() {
            delete innerMatcher;
        }
    public:
        InDelMatcher() {
            innerMatcher = new Matcher();
        }

        Match getMatch(const nrp::NRP *nrp, const nrpsprediction::NRPsPrediction *prediction,
                       const matcher::Score *score) override;


    };
}


#endif //NRPSMATCHER_INDELMATCHER_H
