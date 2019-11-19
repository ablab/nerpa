//
// Created by olga on 13.11.19.
//

#include <Logger/logger.hpp>
#include "OrderedGenesMatcher.h"


class NRP_iterator_rev : public matcher::OrderedGenesMatcher::NRP_iterator {
public:
    explicit NRP_iterator_rev(std::shared_ptr<nrp::NRP> nrp):NRP_iterator(nrp){}

    int getID(int i) override {
        return nrp_->getLen() - i - 1;
    }
};

class NRP_iterator_cycle : public matcher::OrderedGenesMatcher::NRP_iterator {
private:
    int start_ = 0;
public:
    explicit NRP_iterator_cycle(std::shared_ptr<nrp::NRP> nrp, int start):NRP_iterator(nrp){
        start_ = start;
    }

    int getID(int i) override {
        return (start_ + i) % getLen();
    }
};


class NRP_iterator_cycle_rev : public matcher::OrderedGenesMatcher::NRP_iterator {
private:
    int start_ = 0;
public:
    explicit NRP_iterator_cycle_rev(std::shared_ptr<nrp::NRP> nrp, int start):NRP_iterator(nrp){
        start_ = start;
    }

    int getID(int i) override {
        return (start_ - i + getLen()) % getLen();
    }
};


matcher::MatcherBase::Match matcher::OrderedGenesMatcher::getMatch(const std::shared_ptr<nrp::NRP> nrp,
                                                                   const nrpsprediction::NRPsPrediction *prediction,
                                                                   const matcher::Score *score) {
    if (nrp->getType() == nrp::NRP::line) {
        return getLineMatch(false, false, nrp, prediction, score);
    } else if (nrp->getType() == nrp::NRP::cycle) {
        return getCycleMatch(nrp, prediction, score);
    } else {
        return getBranchMatch(nrp, prediction, score);
    }
}


matcher::MatcherBase::Match matcher::OrderedGenesMatcher::getLineMatch(bool can_skip_first, bool can_skip_last,
                                                                       const std::shared_ptr<nrp::NRP> nrp,
                                                                       const nrpsprediction::NRPsPrediction *prediction,
                                                                       const matcher::Score *score) const {
    NRP_iterator* nrp_iterator = new NRP_iterator(nrp);
    matcher::MatcherBase::Match matche1 = getSimpleMatch(can_skip_first, can_skip_last, nrp_iterator, prediction, score);

    NRP_iterator* nrp_iterator_rev = new NRP_iterator_rev(nrp);
    matcher::MatcherBase::Match matche2 = getSimpleMatch(can_skip_first, can_skip_last, nrp_iterator_rev, prediction, score);

    delete nrp_iterator;
    delete nrp_iterator_rev;

    if (matche1.score() > matche2.score()) {
        return matche1;
    }
    return matche2;
}

matcher::MatcherBase::Match matcher::OrderedGenesMatcher::getCycleMatch(const std::shared_ptr<nrp::NRP> nrp,
                                                                        const nrpsprediction::NRPsPrediction *prediction,
                                                                        const matcher::Score *score) const {
    matcher::MatcherBase::Match resMatche(nrp, prediction->getNrpsParts(), score->minScore(nrp->getLen()), score);

    for (int i = 0; i < nrp->getLen(); ++i) {
        NRP_iterator* nrp_iterator = new NRP_iterator_cycle(nrp, i);
        NRP_iterator* nrp_iterator_rev = new NRP_iterator_cycle_rev(nrp, i);

        matcher::MatcherBase::Match matche1 = getSimpleMatch(false, false, nrp_iterator, prediction, score);
        matcher::MatcherBase::Match matche2 = getSimpleMatch(false, false, nrp_iterator_rev, prediction, score);

        delete nrp_iterator;
        delete nrp_iterator_rev;

        if (matche1.score() > resMatche.score()) {
            resMatche = matche1;
        }

        if (matche2.score() > resMatche.score()) {
            resMatche = matche2;
        }
    }

    return resMatche;
}

matcher::MatcherBase::Match matcher::OrderedGenesMatcher::getBranchMatch(const std::shared_ptr<nrp::NRP> nrp,
                                                                         const nrpsprediction::NRPsPrediction *prediction,
                                                                         const matcher::Score *score) const {

    auto matche1 = getLineMatch(true, false, nrp->getLines()[0], prediction, score);
    auto matche2 = getLineMatch(true, false, nrp->getLines()[1], prediction, score);

    if (matche1.score() > matche2.score()) {
        return matche1;
    }
    return matche2;
}

matcher::MatcherBase::Match matcher::OrderedGenesMatcher::getSimpleMatch(bool can_skip_first, bool can_skip_last,
                                                                         matcher::OrderedGenesMatcher::NRP_iterator *nrp_iterator,
                                                                         const nrpsprediction::NRPsPrediction *prediction,
                                                                         const matcher::Score *score) const {
    const auto &nrp_parts = prediction->getNrpsParts();
    std::vector<int> part_id;
    std::vector<int> pos_id;
    size_t plen = 0;

    for (int pid = 0; pid < nrp_parts.size(); ++pid) {
        size_t cur_len = nrp_parts[pid].getAminoacidsPrediction().size();
        for (int i = 0; i < cur_len; ++i) {
            part_id.push_back(pid);
            pos_id.push_back(i);
        }

        plen += cur_len;
    }

    auto get_aapred = [&nrp_parts, &part_id, &pos_id](int ppos) {
        return nrp_parts[part_id[ppos]].getAminoacidsPrediction()[pos_id[ppos]];
    };

    auto nrplen = size_t(nrp_iterator->getLen());

    //dp[nrp_pos][pred_pos] = score for prefix
    std::vector< std::vector<double> > dp(nrplen + 1, std::vector<double>(plen + 1,  score->minScore(nrplen)));
    std::vector< std::vector<int> > px(nrplen + 1, std::vector<int>(plen + 1,  0));
    std::vector< std::vector<int> > py(nrplen + 1, std::vector<int>(plen + 1,  0));

    dp[0][0] = 0;
    for (int npos = 1; npos <= nrplen; ++npos) {
        dp[npos][0] = dp[npos - 1][0] + score->InDelScore();
        px[npos][0] = npos - 1;
        py[npos][0] = 0;
    }

    for (int ppos = 1; ppos <= plen; ++ppos) {
        dp[0][ppos] = dp[0][ppos - 1] + score->InDelScore();
        px[0][ppos] = 0;
        py[0][ppos] = ppos - 1;

        if (dp[0][ppos - pos_id[ppos - 1] - 1] > dp[0][ppos]) {
            dp[0][ppos] = dp[0][ppos - pos_id[ppos - 1] - 1];
            px[0][ppos] = 0;
            py[0][ppos] = ppos - pos_id[ppos - 1] - 1;
        }
    }

    for (int npos = 1; npos <= nrplen; ++npos) {
        for (int ppos = 1; ppos <= plen; ++ppos) {
            dp[npos][ppos] = dp[npos - 1][ppos - 1] + score->aaScore(get_aapred(ppos - 1), nrp_iterator->getAA(npos - 1));
            px[npos][ppos] = npos - 1;
            py[npos][ppos] = ppos - 1;

            if (dp[npos][ppos - pos_id[ppos - 1] - 1] + score->SkipSegment() > dp[npos][ppos]) {
                dp[npos][ppos] = dp[npos][ppos - pos_id[ppos - 1] - 1] + score->SkipSegment();
                px[npos][ppos] = npos;
                py[npos][ppos] = ppos - pos_id[ppos - 1] - 1;
            }

            if (dp[npos][ppos - 1] + score->InDelScore() > dp[npos][ppos]) {
                dp[npos][ppos] = dp[npos][ppos - 1] + score->InDelScore();
                px[npos][ppos] = npos;
                py[npos][ppos] = ppos - 1;
            }

            if (dp[npos - 1][ppos] + score->InDelScore() > dp[npos][ppos]) {
                dp[npos][ppos] = std::max(dp[npos - 1][ppos] + score->InDelScore(), dp[npos][ppos]);
                px[npos][ppos] = npos - 1;
                py[npos][ppos] = ppos;
            }

        }
    }

    matcher::MatcherBase::Match nrPsMatch(nrp_iterator->getNRP(), nrp_parts, dp[nrplen][plen], score);

    int curx = int(nrplen), cury = int(plen);
    while (curx != 0 || cury != 0) {
        if (px[curx][cury] == curx - 1 &&
            py[curx][cury] == cury - 1) {
            nrPsMatch.match(nrp_iterator->getID(curx - 1), part_id[cury - 1], pos_id[cury - 1]);
        }

        int ocurx = curx;
        int ocury = cury;
        curx = px[ocurx][ocury];
        cury = py[ocurx][ocury];
    }

    return nrPsMatch;
}
