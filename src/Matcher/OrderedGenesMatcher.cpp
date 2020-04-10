//
// Created by olga on 13.11.19.
//

#include <Logger/logger.hpp>
#include "OrderedGenesMatcher.h"


class NRP_iterator_skip_first : public matcher::OrderedGenesMatcher::NRP_iterator {
public:
    explicit NRP_iterator_skip_first(std::shared_ptr<nrp::NRP> nrp) : NRP_iterator(nrp) {}

    int getID(int i) override {
        return i + 1;
    }

    int getLen() override {
        return nrp_->getLen() - 1;
    }
};

class NRP_iterator_skip_last : public matcher::OrderedGenesMatcher::NRP_iterator {
public:
    explicit NRP_iterator_skip_last(std::shared_ptr<nrp::NRP> nrp) : NRP_iterator(nrp) {}

    int getLen() override {
        return nrp_->getLen() - 1;
    }
};

class NRP_iterator_skip_tails : public matcher::OrderedGenesMatcher::NRP_iterator {
public:
    explicit NRP_iterator_skip_tails(std::shared_ptr<nrp::NRP> nrp) : NRP_iterator(nrp) {}

    int getID(int i) override {
        return i + 1;
    }

    int getLen() override {
        return nrp_->getLen() - 2;
    }
};

class NRP_iterator_rev : public matcher::OrderedGenesMatcher::NRP_iterator {
public:
    explicit NRP_iterator_rev(std::shared_ptr<nrp::NRP> nrp):NRP_iterator(nrp){}

    int getID(int i) override {
        return nrp_->getLen() - i - 1;
    }
};


class NRP_iterator_rev_skip_first : public matcher::OrderedGenesMatcher::NRP_iterator {
public:
    explicit NRP_iterator_rev_skip_first(std::shared_ptr<nrp::NRP> nrp) : NRP_iterator(nrp) {}

    int getID(int i) override {
        return nrp_->getLen() - i - 1;
    }

    int getLen() override {
        return nrp_->getLen() - 1;
    }
};

class NRP_iterator_rev_skip_last : public matcher::OrderedGenesMatcher::NRP_iterator {
public:
    explicit NRP_iterator_rev_skip_last(std::shared_ptr<nrp::NRP> nrp) : NRP_iterator(nrp) {}

    int getID(int i) override {
        return nrp_->getLen() - i - 2;
    }

    int getLen() override {
        return nrp_->getLen() - 1;
    }
};

class NRP_iterator_rev_skip_tails : public matcher::OrderedGenesMatcher::NRP_iterator {
public:
    explicit NRP_iterator_rev_skip_tails(std::shared_ptr<nrp::NRP> nrp) : NRP_iterator(nrp) {}

    int getID(int i) override {
        return nrp_->getLen() - i - 2;
    }

    int getLen() override {
        return nrp_->getLen() - 2;
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
                                                                   const nrpsprediction::BgcPrediction *prediction,
                                                                   const matcher::Score *score) {

    if (nrp->getType() == nrp::NRP::line) {
        return getLineMatch(true, true, nrp, prediction, score);
    } else if (nrp->getType() == nrp::NRP::cycle) {
        return getCycleMatch(nrp, prediction, score);
    } else {
        return getBranchMatch(nrp, prediction, score);
    }
}


matcher::MatcherBase::Match matcher::OrderedGenesMatcher::getLineMatch(bool can_skip_first, bool can_skip_last,
                                                                       const std::shared_ptr<nrp::NRP> nrp,
                                                                       const nrpsprediction::BgcPrediction *prediction,
                                                                       const matcher::Score *score) const {
    NRP_iterator* nrp_iterator;
    NRP_iterator* nrp_iterator_rev;
    if (can_skip_first && can_skip_last &&
    !nrp->getAminoacid(0).is_AA() && !nrp->getAminoacid(nrp->getLen() - 1).is_AA()) {
        nrp_iterator = new NRP_iterator_skip_tails(nrp);
        nrp_iterator_rev = new NRP_iterator_rev_skip_tails(nrp);
    } else if (can_skip_first && !nrp->getAminoacid(0).is_AA()) {
        nrp_iterator = new NRP_iterator_skip_first(nrp);
        nrp_iterator_rev = new NRP_iterator_rev_skip_first(nrp);
    } else if (can_skip_last && !nrp->getAminoacid(nrp->getLen() - 1).is_AA()) {
        nrp_iterator = new NRP_iterator_skip_last(nrp);
        nrp_iterator_rev = new NRP_iterator_rev_skip_last(nrp);
    } else {
        nrp_iterator = new NRP_iterator(nrp);
        nrp_iterator_rev = new NRP_iterator_rev(nrp);
    }

    matcher::MatcherBase::Match matche1 = getSimpleMatch(nrp_iterator, prediction, score);
    matcher::MatcherBase::Match matche2 = getSimpleMatch(nrp_iterator_rev, prediction, score);

    delete nrp_iterator;
    delete nrp_iterator_rev;

    if (matche1.score() > matche2.score()) {
        return matche1;
    }
    return matche2;
}

matcher::MatcherBase::Match matcher::OrderedGenesMatcher::getCycleMatch(const std::shared_ptr<nrp::NRP> nrp,
                                                                        const nrpsprediction::BgcPrediction *prediction,
                                                                        const matcher::Score *score) const {
    matcher::MatcherBase::Match resMatche(nrp, prediction->getOrfs(), score->minScore(nrp->getLen()), score);

    for (int i = 0; i < nrp->getLen(); ++i) {
        NRP_iterator* nrp_iterator = new NRP_iterator_cycle(nrp, i);
        NRP_iterator* nrp_iterator_rev = new NRP_iterator_cycle_rev(nrp, i);

        matcher::MatcherBase::Match matche1 = getSimpleMatch(nrp_iterator, prediction, score);
        matcher::MatcherBase::Match matche2 = getSimpleMatch(nrp_iterator_rev, prediction, score);

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
                                                                         const nrpsprediction::BgcPrediction *prediction,
                                                                         const matcher::Score *score) const {

    auto matche1 = getLineMatch(true, false, nrp->getLines()[0], prediction, score);
    auto matche2 = getLineMatch(true, false, nrp->getLines()[1], prediction, score);

    if (matche1.score() > matche2.score()) {
        return matche1;
    }
    return matche2;
}

matcher::MatcherBase::Match matcher::OrderedGenesMatcher::getSimpleMatch(matcher::OrderedGenesMatcher::NRP_iterator *nrp_iterator,
                                                                         const nrpsprediction::BgcPrediction *prediction,
                                                                         const matcher::Score *score) const {
    const auto &nrp_parts = prediction->getOrfs();
    std::vector<int> part_id;
    std::vector<int> pos_id;
    size_t plen = 0;

    for (int pid = 0; pid < nrp_parts.size(); ++pid) {
        size_t cur_len = nrp_parts[pid].getAAdomainPrediction().size();
        for (int i = 0; i < cur_len; ++i) {
            part_id.push_back(pid);
            pos_id.push_back(i);
        }

        plen += cur_len;
    }

    auto get_aapred = [&nrp_parts, &part_id, &pos_id](int ppos) {
        return nrp_parts[part_id[ppos]].getAAdomainPrediction()[pos_id[ppos]];
    };

    auto nrplen = size_t(nrp_iterator->getLen());

    //dp[nrp_pos][pred_pos] = score for prefix
    std::vector<std::vector< std::vector<double> > > dp(2,
                                                        std::vector<std::vector<double>>(nrplen + 1,
                                                        std::vector<double>(plen + 1,  score->minScore(nrplen))));
    std::vector<std::vector< std::vector<int> > > pgap(2,
                                                       std::vector<std::vector<int> >(nrplen + 1,
                                                       std::vector<int>(plen + 1,  0)));
    std::vector<std::vector< std::vector<int> > > px(2, std::vector<std::vector<int> >(nrplen + 1,
                                                        std::vector<int>(plen + 1,  0)));
    std::vector<std::vector< std::vector<int> > > py(2, std::vector<std::vector<int> >(nrplen + 1,
                                                        std::vector<int>(plen + 1,  0)));

    dp[0][0][0] = 0;
    dp[1][0][0] = -1;
    for (int npos = 1; npos <= nrplen; ++npos) {
        dp[0][npos][0] = dp[0][npos - 1][0] + score->InsertionScore();
        pgap[0][npos][0] = 0;
        px[0][npos][0] = npos - 1;
        py[0][npos][0] = 0;

        if (dp[1][npos - 1][0] + score->InsertionScore() > dp[0][npos][0]) {
            pgap[0][npos][0] = 1;
        }

        dp[1][npos][0] = dp[1][npos - 1][0] + score->continueGap();
        pgap[1][npos][0] = 1;
        px[1][npos][0] = npos - 1;
        py[1][npos][0] = 0;

        if (dp[0][npos - 1][0] + score->openGap() > dp[1][npos][0]) {
            pgap[1][npos][0] = 0;
        }
    }

    for (int ppos = 1; ppos <= plen; ++ppos) {
        dp[0][0][ppos] = dp[0][0][ppos - 1] + score->DeletionScore();
        pgap[0][0][ppos] = 0;
        px[0][0][ppos] = 0;
        py[0][0][ppos] = ppos - 1;

        if ((ppos == plen || part_id[ppos] != part_id[ppos - 1]) &&
        dp[0][0][ppos - pos_id[ppos - 1] - 1] + score->SkipSegment() > dp[0][0][ppos]) {
            dp[0][0][ppos] = dp[0][0][ppos - pos_id[ppos - 1] - 1] + score->SkipSegment();
            pgap[0][0][ppos] = 0;
            px[0][0][ppos] = 0;
            py[0][0][ppos] = ppos - pos_id[ppos - 1] - 1;
        }
    }

    for (int npos = 1; npos <= nrplen; ++npos) {
        for (int ppos = 1; ppos <= plen; ++ppos) {
            dp[0][npos][ppos] = dp[0][npos - 1][ppos - 1] + score->aaScore(get_aapred(ppos - 1), nrp_iterator->getAA(npos - 1));
            pgap[0][npos][ppos] = 0;
            px[0][npos][ppos] = npos - 1;
            py[0][npos][ppos] = ppos - 1;

            if ((ppos == plen || part_id[ppos] != part_id[ppos - 1]) &&
            dp[0][npos][ppos - pos_id[ppos - 1] - 1] + score->SkipSegment() > dp[0][npos][ppos]) {
                dp[0][npos][ppos] = dp[0][npos][ppos - pos_id[ppos - 1] - 1] + score->SkipSegment();
                pgap[0][npos][ppos] = 0;
                px[0][npos][ppos] = npos;
                py[0][npos][ppos] = ppos - pos_id[ppos - 1] - 1;
            }

            if (dp[0][npos][ppos - 1] + score->DeletionScore() > dp[0][npos][ppos]) {
                dp[0][npos][ppos] = dp[0][npos][ppos - 1] + score->DeletionScore();
                pgap[0][npos][ppos] = 0;
                px[0][npos][ppos] = npos;
                py[0][npos][ppos] = ppos - 1;
            }

            if (dp[0][npos - 1][ppos] + score->InsertionScore() > dp[0][npos][ppos]) {
                dp[0][npos][ppos] = std::max(dp[0][npos - 1][ppos] + score->InsertionScore(), dp[0][npos][ppos]);
                pgap[0][npos][ppos] = 0;
                px[0][npos][ppos] = npos - 1;
                py[0][npos][ppos] = ppos;
            }

            if (ppos == plen || pos_id[ppos] == 0) {
                dp[1][npos][ppos] = dp[1][npos - 1][ppos] + score->continueGap();
                pgap[1][npos][ppos] = 1;
                px[1][npos][ppos] = npos - 1;
                py[1][npos][ppos] = ppos;

                if (dp[0][npos - 1][ppos] + score->openGap() > dp[1][npos][ppos]) {
                    dp[1][npos][ppos] = dp[0][npos - 1][ppos] + score->openGap();
                    pgap[1][npos][ppos] = 0;
                }

                if (dp[0][npos][ppos] < dp[1][npos][ppos]) {
                    dp[0][npos][ppos] = dp[1][npos][ppos];
                    pgap[0][npos][ppos] = 1;
                    px[0][npos][ppos] = npos;
                    py[0][npos][ppos] = ppos;
                }
            }

        }
    }

    matcher::MatcherBase::Match nrPsMatch(nrp_iterator->getNRP(), nrp_parts, dp[0][nrplen][plen], score);

    int curx = int(nrplen), cury = int(plen), curgap = 0;
    while (curx != 0 || cury != 0) {
        if (px[curgap][curx][cury] == curx - 1 &&
            py[curgap][curx][cury] == cury - 1) {
            nrPsMatch.match(nrp_iterator->getID(curx - 1), part_id[cury - 1], pos_id[cury - 1]);
        }

        int ocurx = curx;
        int ocury = cury;
        curx = px[curgap][ocurx][ocury];
        cury = py[curgap][ocurx][ocury];
        curgap = pgap[curgap][ocurx][ocury];
    }

    nrPsMatch.setScore(score->resultScore(dp[0][nrplen][plen], nrplen));
    return nrPsMatch;
}
