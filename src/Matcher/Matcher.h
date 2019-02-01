//
// Created by olga on 22.01.19.
//

#ifndef NRPSMATCHER_MATCHER_H
#define NRPSMATCHER_MATCHER_H

#include <NRPsPrediction/AminoacidPrediction.h>
#include <NRP/NRP.h>
#include <Matcher/Score/Score.h>
#include "Segment.h"

namespace matcher {
    typedef aminoacid::Aminoacid aacid;
    class Matcher {
    private:
        const nrp::NRP& nrp;
        const nrpsprediction::NRPsPrediction& prediction;
        const Score score = Score();
    public:
        class Match {
        private:
            const nrp::NRP* nrp;
            double scr;
            std::vector<nrpsprediction::NRPsPart> nrpParts;
            std::vector<int> parts_id;
            std::vector<int> parts_pos;
        public:
            Match() = default;

            Match(const nrp::NRP* nrp, std::vector<nrpsprediction::NRPsPart> nrpParts, double scr):
                    nrp(nrp), nrpParts(nrpParts), scr(scr) {
                parts_id.resize(nrp->getLen(), -1);
                parts_pos.resize(nrp->getLen(), -1);
            }

            void match(int pos, int part_id, int part_pos);
            void print(std::ofstream& out);
            void print_short(std::ofstream& out);
            void print_short_prediction(std::ofstream& out);
            void print_csv(std::ofstream& out);
            double score();
            bool isMatched(int i);
            std::vector<std::pair<int, int> > getMatchs();

            bool operator < (Match b);
        };
        Matcher(const nrp::NRP &nrp, const nrpsprediction::NRPsPrediction& prediction):
                nrp(nrp), prediction(prediction) {}

        matcher::Matcher::Match getMatch() const;
        std::vector<Segment> matche_seg(const nrpsprediction::NRPsPart& predict_part) const;
    private:
        matcher::Matcher::Match getLineMatch(bool can_skip_first = true, bool can_skip_last = true) const;
        matcher::Matcher::Match getCycleMatch() const;
        matcher::Matcher::Match getBranchMatch() const;


        matcher::Matcher::Match updateMatch(const nrpsprediction::NRPsPrediction& nrPsPrediction, matcher::Matcher::Match match, int bg) const;
        matcher::Matcher::Match isCoverLine(std::vector<Segment>& segments,
                                  const std::vector<int>& toSmallId, const std::vector<int>& toBigId, int len) const;

        std::vector<aacid> getSubset(std::vector<aacid> vector, int l, int r, int stp) const;

        int addSegments(const std::vector<Segment> &vector, std::vector<Segment> &segments, bool skip_first, bool skip_last, int id) const;
    };
}


#endif //NRPSMATCHER_MATCHER_H
