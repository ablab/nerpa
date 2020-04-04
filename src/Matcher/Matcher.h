//
// Created by olga on 22.01.19.
//

#ifndef NRPSMATCHER_MATCHER_H
#define NRPSMATCHER_MATCHER_H

#include <NRPsPrediction/AAdomain_Prediction.h>
#include <NRP/NRP.h>
#include <Matcher/Score/Base/Score.h>

#include <utility>
#include "Segment.h"
#include "MatcherBase.h"

namespace matcher {
    typedef aminoacid::Aminoacid aacid;
    class Matcher : public MatcherBase {
    protected:
        std::shared_ptr<nrp::NRP> nrp;
        const nrpsprediction::BGC_Prediction* prediction{};
        const Score* score{};
    public:
        Matcher(std::shared_ptr<nrp::NRP> nrp,
                const nrpsprediction::BGC_Prediction* prediction,
                const Score* score):
                nrp(nrp), prediction(prediction), score(score) {
        }

        Matcher() = default;

        matcher::MatcherBase::Match getMatch() const;
        matcher::MatcherBase::Match getMatch(std::shared_ptr<nrp::NRP> nrp,
                                             const nrpsprediction::BGC_Prediction* prediction,
                                             const Score* score) override;
        std::vector<Segment> matche_seg(const int part_id) const;

    protected:
        matcher::MatcherBase::Match getLineMatch(bool can_skip_first = true, bool can_skip_last = true) const;
        matcher::MatcherBase::Match getCycleMatch() const;
        matcher::MatcherBase::Match getBranchMatch() const;

        virtual matcher::MatcherBase::Match updateMatch(const nrpsprediction::BGC_Prediction& nrPsPrediction,
                                            matcher::MatcherBase::Match match, int bg,
                                            std::vector<Segment>& matched_parts_id) const;
        matcher::MatcherBase::Match setUpdateMatch(const nrpsprediction::BGC_Prediction& nrPsPrediction,
                                                matcher::MatcherBase::Match match, int bg,
                                                std::vector<Segment>& matched_parts_id) const;
        matcher::MatcherBase::Match isCoverLine(std::vector<Segment>& segments,
                                  const std::vector<int>& toSmallId, const std::vector<int>& toBigId, int len,
                                            std::vector<Segment>& matched_parts_id) const;

        std::vector<aacid> getSubset(std::vector<aacid> vector, int l, int r, int stp) const;

        int addSegments(const std::vector<Segment> &vector, std::vector<Segment> &segments, bool skip_first, bool skip_last, int id) const;
    };
}


#endif //NRPSMATCHER_MATCHER_H
