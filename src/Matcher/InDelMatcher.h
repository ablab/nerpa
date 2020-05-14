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
        bool insertion = true;
        bool deletion = true;
        MatcherBase* innerMatcher = nullptr;

        ~InDelMatcher() {
            delete innerMatcher;
        }

        Match getDeleteMatch(std::shared_ptr<nrp::NRP> nrp, const nrpsprediction::BgcPrediction *prediction,
                             const matcher::Score *score);

        Match getInsertMatch(std::shared_ptr<nrp::NRP> nrp, const nrpsprediction::BgcPrediction *prediction,
                             const matcher::Score *score);
    public:
        InDelMatcher() {
            innerMatcher = new Matcher();
        }

        explicit InDelMatcher(MatcherBase* matcher) {
            innerMatcher = matcher;
        }

        InDelMatcher(bool insertion, bool deletion) {
            innerMatcher = new Matcher();
            this->insertion = insertion;
            this->deletion = deletion;
        }


        InDelMatcher(MatcherBase* matcher, bool insertion, bool deletion) {
            innerMatcher = matcher;
            this->insertion = insertion;
            this->deletion = deletion;
        }

        Match getMatch(std::shared_ptr<nrp::NRP> nrp, const nrpsprediction::BgcPrediction *prediction,
                       const matcher::Score *score) override;
    };
}


#endif //NRPSMATCHER_INDELMATCHER_H
