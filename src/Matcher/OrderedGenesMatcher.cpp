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

    auto is_rep_aa = [&nrp_parts, &part_id, &pos_id](int ppos) {
        return nrp_parts[part_id[ppos]].getAAdomainPrediction()[pos_id[ppos]].is_rep();
    };

    auto is_rep_orf = [&nrp_parts, &part_id, &pos_id](int ppos) {
        return nrp_parts[part_id[ppos]].is_rep();
    };

    auto orf_len = [&nrp_parts, &part_id, &pos_id](int ppos) {
        return nrp_parts[part_id[ppos]].getAAdomainPrediction().size();
    };

    auto nrplen = size_t(nrp_iterator->getLen());

    //dp[nrp_pos][pred_pos] = score for prefix
    std::vector< std::vector<double> > dp(std::vector<std::vector<double>>(nrplen + 1,
            std::vector<double>(plen + 1,  score->minScore(nrplen))));
    std::vector< std::vector<int> > px(std::vector<std::vector<int> >(nrplen + 1,
            std::vector<int>(plen + 1,  0)));
    std::vector< std::vector<int> > py(std::vector<std::vector<int> >(nrplen + 1,
            std::vector<int>(plen + 1,  0)));
    std::vector< std::vector<int> > pmtc(std::vector<std::vector<int> >(nrplen + 1,
            std::vector<int>(plen + 1,  0)));

    dp[0][0] = 0;
    for (int npos = 1; npos <= nrplen; ++npos) {
        dp[npos][0] = dp[npos - 1][0] + score->InsertionScore();
        px[npos][0] = npos - 1;
        py[npos][0] = 0;
    }

    for (int ppos = 1; ppos <= plen; ++ppos) {
        dp[0][ppos] = dp[0][ppos - 1] + score->DeletionScore();
        px[0][ppos] = 0;
        py[0][ppos] = ppos - 1;

        if ((ppos == plen || part_id[ppos] != part_id[ppos - 1]) &&
        dp[0][ppos - pos_id[ppos - 1] - 1] + score->SkipSegment() > dp[0][ppos]) {
            dp[0][ppos] = dp[0][ppos - pos_id[ppos - 1] - 1] + score->SkipSegment();
            px[0][ppos] = 0;
            py[0][ppos] = ppos - pos_id[ppos - 1] - 1;
        }
    }

    for (int npos = 1; npos <= nrplen; ++npos) {
        for (int ppos = 1; ppos <= plen; ++ppos) {
            double match_score =  score->aaScore(get_aapred(ppos - 1), nrp_iterator->getAA(npos - 1));
            dp[npos][ppos] = dp[npos - 1][ppos - 1] + match_score;
            px[npos][ppos] = npos - 1;
            py[npos][ppos] = ppos - 1;
            pmtc[npos][ppos] = 1;

            if (is_rep_aa(ppos - 1) &&
                dp[npos - 1][ppos] + match_score > dp[npos][ppos]) {
                dp[npos][ppos] = dp[npos - 1][ppos] + match_score;
                px[npos][ppos] = npos - 1;
                py[npos][ppos] = ppos;
                pmtc[npos][ppos] = 1;
            }

            if (is_rep_orf(ppos - 1) && pos_id[ppos - 1] == 0 &&
            dp[npos - 1][ppos + orf_len(ppos - 1) - 1] + match_score > dp[npos][ppos]
            ) {
                dp[npos][ppos] = dp[npos - 1][ppos + orf_len(ppos - 1) - 1] + match_score;
                px[npos][ppos] = npos - 1;
                py[npos][ppos] = ppos + orf_len(ppos - 1) - 1 ;
                pmtc[npos][ppos] = 1;
            }

            if ((ppos == plen || part_id[ppos] != part_id[ppos - 1]) &&
            dp[npos][ppos - pos_id[ppos - 1] - 1] + score->SkipSegment() > dp[npos][ppos]) {
                dp[npos][ppos] = dp[npos][ppos - pos_id[ppos - 1] - 1] + score->SkipSegment();
                px[npos][ppos] = npos;
                py[npos][ppos] = ppos - pos_id[ppos - 1] - 1;
                pmtc[npos][ppos] = 0;
            }

            if (dp[npos][ppos - 1] + score->DeletionScore() > dp[npos][ppos]) {
                dp[npos][ppos] = dp[npos][ppos - 1] + score->DeletionScore();
                px[npos][ppos] = npos;
                py[npos][ppos] = ppos - 1;
                pmtc[npos][ppos] = 0;
            }

            if (dp[npos - 1][ppos] + score->InsertionScore() > dp[npos][ppos]) {
                dp[npos][ppos] = std::max(dp[npos - 1][ppos] + score->InsertionScore(), dp[npos][ppos]);
                px[npos][ppos] = npos - 1;
                py[npos][ppos] = ppos;
                pmtc[npos][ppos] = 0;
            }
        }
    }

    matcher::MatcherBase::Match nrPsMatch(nrp_iterator->getNRP(), nrp_parts, dp[nrplen][plen], score);

    int curx = int(nrplen), cury = int(plen);
    while (curx != 0 || cury != 0) {
        if (pmtc[curx][cury]) {
            nrPsMatch.match(nrp_iterator->getID(curx - 1), part_id[cury - 1], pos_id[cury - 1]);
        }

        int ocurx = curx;
        int ocury = cury;
        curx = px[ocurx][ocury];
        cury = py[ocurx][ocury];
    }

    nrPsMatch.setScore(score->resultScore(dp[nrplen][plen], nrplen));
    return nrPsMatch;
}