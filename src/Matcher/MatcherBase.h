//
// Created by olga on 03.07.19.
//

#ifndef NRPSMATCHER_MATCHERBASE_H
#define NRPSMATCHER_MATCHERBASE_H

#include <NRPsPrediction/AminoacidPrediction.h>
#include <NRP/NRP.h>
#include <Matcher/Score/Base/Score.h>

#include <utility>
#include "Segment.h"

namespace matcher {
    typedef aminoacid::Aminoacid aacid;
    class MatcherBase {
    public:
        class Match {
        private:
            std::shared_ptr<nrp::NRP> nrp;
            const Score* scoreFun;
            double scr;
            std::vector<nrpsprediction::NRPsPart> nrpParts;
            std::vector<int> parts_id;
            std::vector<int> parts_pos;
        public:
            Match(std::shared_ptr<nrp::NRP> nrp, std::vector<nrpsprediction::NRPsPart> nrpParts, double scr, const Score* score):
                    nrp(nrp), nrpParts(std::move(nrpParts)), scr(scr), scoreFun(score) {
                parts_id.resize(nrp->getFullLen(), -1);
                parts_pos.resize(nrp->getFullLen(), -1);
            }

            void match(int pos, int part_id, int part_pos);
            void print(std::ofstream& out);
            void print_short(std::ofstream& out);
            void print_short_prediction(std::ofstream& out);
            void print_csv(std::ofstream& out);
            double score() const;
            void setScore(double score);
            int getCntMatch();
            bool isMatched(int i);
            std::vector<std::pair<int, int> > getMatchs();

            bool operator < (Match b) const;
        };

        MatcherBase() = default;

        virtual matcher::MatcherBase::Match getMatch(const std::shared_ptr<nrp::NRP> nrp,
                                                     const nrpsprediction::NRPsPrediction *prediction,
                                                     const matcher::Score *score) = 0;
    };
}

#endif //NRPSMATCHER_MATCHERBASE_H
