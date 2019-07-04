//
// Created by olga on 22.01.19.
//

#ifndef NRPSMATCHER_MATCHER_H
#define NRPSMATCHER_MATCHER_H

#include <NRPsPrediction/AminoacidPrediction.h>
#include <NRP/NRP.h>
#include <Matcher/Score/Base/Score.h>

#include <utility>
#include "Segment.h"
#include "MatcherBase.h"

namespace matcher {
    typedef aminoacid::Aminoacid aacid;
    class Matcher : public MatcherBase {
    private:
        const nrp::NRP& nrp;
        const nrpsprediction::NRPsPrediction& prediction;
        const Score* score;
    public:
        Matcher(const nrp::NRP &nrp, const nrpsprediction::NRPsPrediction& prediction, const Score* score):
                MatcherBase(nrp, prediction, score), nrp(nrp), prediction(prediction), score(score) {
        }

        matcher::MatcherBase::Match getMatch() const;
        std::vector<Segment> matche_seg(const int part_id) const;
    private:
        matcher::MatcherBase::Match getLineMatch(bool can_skip_first = true, bool can_skip_last = true) const;
        matcher::MatcherBase::Match getCycleMatch() const;
        matcher::MatcherBase::Match getBranchMatch() const;


        matcher::MatcherBase::Match updateMatch(const nrpsprediction::NRPsPrediction& nrPsPrediction,
                                            matcher::MatcherBase::Match match, int bg,
                                            std::vector<Segment>& matched_parts_id) const;
        matcher::MatcherBase::Match isCoverLine(std::vector<Segment>& segments,
                                  const std::vector<int>& toSmallId, const std::vector<int>& toBigId, int len,
                                            std::vector<Segment>& matched_parts_id) const;

        std::vector<aacid> getSubset(std::vector<aacid> vector, int l, int r, int stp) const;

        int addSegments(const std::vector<Segment> &vector, std::vector<Segment> &segments, bool skip_first, bool skip_last, int id) const;
        void matchSingleUnits(Match& match, std::vector<bool>& used_pos) const;
    };
}


#endif //NRPSMATCHER_MATCHER_H
