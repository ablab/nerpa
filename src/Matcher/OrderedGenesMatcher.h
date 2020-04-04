//
// Created by olga on 13.11.19.
//

#ifndef NERPA_ORDEREDGENESMATCHER_H
#define NERPA_ORDEREDGENESMATCHER_H


#include "MatcherBase.h"

namespace matcher {
    class OrderedGenesMatcher : public MatcherBase {
    public:
        class NRP_iterator {
        protected:
            std::shared_ptr<nrp::NRP> nrp_;
        public:
            explicit NRP_iterator(std::shared_ptr<nrp::NRP> nrp): nrp_(nrp){}
            aminoacid::Aminoacid getAA(int i) {
                return nrp_->getAminoacid(getID(i));
            }

            virtual int getID(int i) {
                return i;
            }

            int getLen() {
                return nrp_->getLen();
            }

            std::shared_ptr<nrp::NRP> getNRP() {
                return nrp_;
            }
        };

    private:
        matcher::MatcherBase::Match getSimpleMatch(bool can_skip_first, bool can_skip_last,
                                                   NRP_iterator* nrp_iterator,
                                                   const nrpsprediction::BgcPrediction *prediction,
                                                   const matcher::Score *score) const;

        matcher::MatcherBase::Match getLineMatch(bool can_skip_first, bool can_skip_last,
                                                 std::shared_ptr<nrp::NRP> nrp, const nrpsprediction::BgcPrediction *prediction,
                                                 const matcher::Score *score) const;
        matcher::MatcherBase::Match getCycleMatch(std::shared_ptr<nrp::NRP> nrp, const nrpsprediction::BgcPrediction *prediction,
                                                  const matcher::Score *score) const;
        matcher::MatcherBase::Match getBranchMatch(std::shared_ptr<nrp::NRP> nrp, const nrpsprediction::BgcPrediction *prediction,
                                                   const matcher::Score *score) const;
    public:
        Match getMatch(std::shared_ptr<nrp::NRP> nrp, const nrpsprediction::BgcPrediction *prediction,
                       const matcher::Score *score) override;

    };
}


#endif //NERPA_ORDEREDGENESMATCHER_H
