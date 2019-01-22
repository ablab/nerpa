//
// Created by olga on 22.01.19.
//

#ifndef NRPSMATCHER_MATCHER_H
#define NRPSMATCHER_MATCHER_H

#include <NRPsPrediction/AminoacidPrediction.h>
#include <NRP/NRP.h>
#include <Matcher/Score/Score.h>

namespace matcher {
    class Matcher {
    private:
        const nrp::NRP& nrp;
        const nrpsprediction::NRPsPrediction& prediction;
        const Score score = Score();
    public:
        class Match {
        private:
            const nrp::NRP* nrp;
            std::vector<nrpsprediction::NRPsPart> nrpParts;
            std::vector<int> parts_id;
            std::vector<int> parts_pos;
        public:
            Match() = default;

            Match(const nrp::NRP* nrp, std::vector<nrpsprediction::NRPsPart> nrpParts): nrp(nrp), nrpParts(nrpParts) {
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
    private:
        matcher::Matcher::Match getLineMatch() const;
        matcher::Matcher::Match getCycleMatch() const;
        matcher::Matcher::Match getBranchMatch() const;


        matcher::Matcher::Match updateMatch(const nrpsprediction::NRPsPrediction& nrPsPrediction, matcher::Matcher::Match match, int bg) const;
        virtual matcher::Matcher::Match isCoverLine(std::vector<nrp::NRP::Segment>& segments,
                                  const std::vector<int>& toSmallId, const std::vector<int>& toBigId) const;
    };
}


#endif //NRPSMATCHER_MATCHER_H
