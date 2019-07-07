//
// Created by olga on 22.01.19.
//

#include <NRP/NRPLine.h>
#include <algorithm>
#include <fstream>
#include "Matcher.h"

matcher::MatcherBase::Match matcher::Matcher::getMatch() const {
    if (nrp->getType() == nrp::NRP::line) {
        return getLineMatch(true, true);
    } else if (nrp->getType() == nrp::NRP::cycle) {
        return getCycleMatch();
    } else {
        return getBranchMatch();
    }
}

matcher::MatcherBase::Match matcher::Matcher::getLineMatch(bool can_skip_first, bool can_skip_last) const {
    std::vector<Segment> segments;
    auto nrpparts = prediction->getNrpsParts();
    std::vector<int> toSmallId(nrpparts.size(), -1);
    std::vector<int> toBigId;
    int len = nrp->getLen();
    bool skip_first = (can_skip_first && len > 0 && (nrp->getAminoacid(0).get_name() == "none"));
    bool skip_last = (can_skip_last && len > 0 && (nrp->getAminoacid(len - 1).get_name() == "none"));

    len -= (skip_first + skip_last);

    for (int i = 0; i < nrpparts.size(); ++i) {
        std::vector<Segment> part_seg = matche_seg(i);

        int cnt_add = addSegments(part_seg, segments, skip_first, skip_last, i);

        if (cnt_add >= 2) {
            toSmallId[i] = toBigId.size();
            toBigId.push_back(i);
        }
    }
    std::vector<Segment> matched_parts_id;
    return updateMatch(*prediction, isCoverLine(segments, toSmallId, toBigId, len, matched_parts_id), skip_first, matched_parts_id);
}

matcher::MatcherBase::Match matcher::Matcher::getCycleMatch() const {
    std::vector<Segment> segments;
    auto nrpparts = prediction->getNrpsParts();
    std::vector<int> toSmallId(nrpparts.size(), -1);
    std::vector<int> toBigId;
    int len = nrp->getLen();

    for (int i = 0; i < nrpparts.size(); ++i) {
        std::vector<Segment> part_seg = matche_seg(i);

        for (int j = 0; j < part_seg.size(); ++j) {
            segments.push_back(Segment(part_seg[j].l, part_seg[j].r, i, part_seg[j].rev, part_seg[j].scor));
        }

        if (part_seg.size() >= 2) {
            toSmallId[i] = toBigId.size();
            toBigId.push_back(i);
        }
    }

    double best_score = score->minScore(len);
    matcher::MatcherBase::Match resMatchs(nrp.get(), nrpparts, score->minScore(len), score);

    for (int bg = 0; bg < len; ++bg) {
        std::vector<Segment> tmpSeg;
        for (int i = 0; i < segments.size(); ++i) {
            tmpSeg.push_back(Segment(segments[i].l - bg, segments[i].r - bg, segments[i].part_id, segments[i].rev, segments[i].scor));
            if (tmpSeg[tmpSeg.size() - 1].l < 0) {
                tmpSeg[tmpSeg.size() - 1].l += len;
                tmpSeg[tmpSeg.size() - 1].r += len;
            }
        }

        std::vector<Segment> matched_parts_id;
        matcher::MatcherBase::Match curMatch = updateMatch(*prediction, isCoverLine(tmpSeg, toSmallId, toBigId, len, matched_parts_id), bg, matched_parts_id);
        if (curMatch.score() > best_score) {
            resMatchs = curMatch;
            best_score = curMatch.score();
        }
    }


    return resMatchs;
}


matcher::MatcherBase::Match matcher::Matcher::getBranchMatch() const {
    std::vector<std::shared_ptr<nrp::NRP>> lines = nrp->getLines();

    assert(lines[0]->getFormula(0) == lines[1]->getFormula(0));

    Matcher matcher1(lines[0], prediction, score);
    Matcher matcher2(lines[1], prediction, score);
    matcher::MatcherBase::Match m1 = matcher1.getLineMatch(true, false);
    matcher::MatcherBase::Match m2 = matcher2.getLineMatch(true, false);
    if (m1.score() > m2.score()) {
        return m1;
    } else {
        return m2;
    }
}

matcher::MatcherBase::Match
matcher::Matcher::updateMatch(const nrpsprediction::NRPsPrediction &nrPsPrediction, matcher::MatcherBase::Match match,
                              int bg, std::vector<Segment>& matched_parts_id) const {
    auto parts = nrPsPrediction.getNrpsParts();
    auto short_parts = nrPsPrediction.getShortParts();
    parts.insert(parts.end(), short_parts.begin(), short_parts.end());
    matcher::MatcherBase::Match nmatch(nrp.get(), parts, match.score(), score);
    std::vector<std::pair<int, int> > part_id_pos = match.getMatchs();
    std::vector<bool> matched_nrp_pos(nrp->getLen(), false);
    for (int i = 0; i < part_id_pos.size(); ++i) {
        nmatch.match((i + bg)%nrp->getLen(), part_id_pos[i].first, part_id_pos[i].second);
        if (part_id_pos[i].first > -1) {
            matched_nrp_pos[(i + bg) % nrp->getLen()] = true;
        }
    }

    for (int i = 0; i < matched_parts_id.size(); ++i) {
        matched_parts_id[i].l += bg;
        matched_parts_id[i].r += bg;
        matched_parts_id[i].l %= nrp->getLen();
        matched_parts_id[i].r %= nrp->getLen();
    }

    std::sort(matched_parts_id.begin(), matched_parts_id.end(), [](Segment a, Segment b) -> bool {
        return a.part_id < b.part_id;
    });
    matched_parts_id.resize(std::unique(matched_parts_id.begin(), matched_parts_id.end(), [](Segment a, Segment b) -> bool {
        return a.part_id == b.part_id;
    }) - matched_parts_id.begin());

    nmatch.setScore(match.score());
    matchSingleUnits(nmatch, matched_nrp_pos);
    nmatch.setScore(score->resultScore(nmatch.score(), nrp->getLen(), matched_parts_id, *prediction, *nrp));
    return nmatch;
}

matcher::MatcherBase::Match
matcher::Matcher::isCoverLine(std::vector<Segment> &segments,
                      const std::vector<int> &toSmallId, const std::vector<int> &toBigId, int len, std::vector<Segment>& matched_parts_id) const {

    if (toBigId.size() > 20) {
        return matcher::MatcherBase::Match(nrp.get(), prediction->getNrpsParts(), -1, score);
    }

    std::sort(segments.begin(), segments.end());

    std::vector<std::vector<std::vector<double> > > d(len + 1,
                                                      std::vector<std::vector<double> >((1 << toBigId.size()),
                                                                                        std::vector<double>(2, score->minScore(len))));
    std::vector<std::vector<std::pair<int, int> > > p(len + 1,
                                                      std::vector<std::pair<int, int> >((1 << toBigId.size()),
                                                                                        std::make_pair(-1, -1)));

    std::vector<std::vector<std::vector<int> > > pgp(len + 1, std::vector<std::vector<int> >((1 << toBigId.size()),
                                                                                             std::vector<int>(2, 0)));
    std::vector<std::vector<int> > pa(len + 1, std::vector<int>((1 << toBigId.size()), -1));
    d[0][0][1] = 0;

    int curseg = 0;

    for (int pos = 0; pos <= len; ++pos) {
        int lstseg = curseg;
        for (int msk = 0; msk < (1 << (toBigId.size())); ++msk) {
            if (pos != 0 && d[pos][msk][0] < d[pos - 1][msk][1] + score->openGap()) {
                d[pos][msk][0] = d[pos - 1][msk][1] + score->openGap();
                pgp[pos][msk][0] = 1;
            }
            if (pos != 0 && d[pos][msk][0] < d[pos - 1][msk][0] + score->continueGap()) {
                d[pos][msk][0] = d[pos - 1][msk][0] + score->continueGap();
                pgp[pos][msk][0] = 0;
            }
            curseg = lstseg;

            while (curseg < segments.size() && segments[curseg].l == pos) {
                if (segments[curseg].r >= len) {
                    ++curseg;
                    continue;
                }
                int nmsk = msk;

                if (toSmallId[segments[curseg].part_id] != -1) {
                    if (((msk >> (toSmallId[segments[curseg].part_id])) & 1) != 0) {
                        ++curseg;
                        continue;
                    }

                    nmsk = msk | (1 << (toSmallId[segments[curseg].part_id]));
                }

                for (int gp = 0; gp < 2; ++gp) {
                    if (d[segments[curseg].r + 1][nmsk][1] < d[pos][msk][gp] + score->addSegment(segments[curseg])) {
                        d[segments[curseg].r + 1][nmsk][1] = d[pos][msk][gp] + score->addSegment(segments[curseg]);
                        p[segments[curseg].r + 1][nmsk].first = pos;
                        p[segments[curseg].r + 1][nmsk].second = msk;
                        pa[segments[curseg].r + 1][nmsk] = curseg;
                        pgp[segments[curseg].r + 1][nmsk][1] = gp;
                    }
                }

                ++curseg;
            }
        }
    }

    double mn = score->minScore(len);
    int rmsk = 0;
    int gp = 0;
    for (int msk = 0; msk < (1 << (toBigId.size())); ++msk) {
        if (mn <= d[len][msk][0]) {
            mn = d[len][msk][0];
            rmsk = msk;
            gp = 0;
        }

        if (mn <= d[len][msk][1]) {
            mn = d[len][msk][1];
            rmsk = msk;
            gp = 1;
        }
    }

    auto parts = prediction->getNrpsParts();
    auto short_parts = prediction->getShortParts();
    parts.insert(parts.end(), short_parts.begin(), short_parts.end());
    matcher::MatcherBase::Match nrPsMatch(nrp.get(), parts, mn, score);
    int pos = len;
    while (pos > 0) {
        int nxtp = p[pos][rmsk].first;
        int nxtmsk = p[pos][rmsk].second;
        int ngp = pgp[pos][rmsk][gp];
        if (gp == 0) {
            nxtp = pos - 1;
            nxtmsk = rmsk;
        }

        int seg = pa[pos][rmsk];
        if (gp == 1 && seg != -1) {
            int j = segments[seg].r - segments[seg].l;
            if (segments[seg].rev) {
                j = 0;
            }

            int curp = pos - 1;
            while (curp >= nxtp) {
                nrPsMatch.match(curp, segments[seg].part_id, j);
                matched_parts_id.push_back(segments[seg]);
                if (segments[seg].rev) {
                    ++j;
                } else {
                    --j;
                }
                --curp;
            }
        }

        pos = nxtp;
        rmsk = nxtmsk;
        gp = ngp;
    }

    return nrPsMatch;
}


std::vector<matcher::Segment> matcher::Matcher::matche_seg(const int part_id) const {
    auto parts = prediction->getNrpsParts();
    nrpsprediction::NRPsPart predict_part = prediction->getNrpsParts()[part_id];
    std::vector<Segment> segs;
    std::vector<aacid> amns = nrp->getAminoacids();
    amns.resize(nrp->getLen());
    int part_len = predict_part.getAminoacidsPrediction().size();

    if (part_len > amns.size() || part_len == 0) {
        return segs;
    }

    for (int l = 0; l < amns.size(); ++l) {
        for (int stp = -1; stp < 2; stp += 2) {
            int r =  (l + stp*(part_len - 1) + amns.size()) % amns.size();
            if (nrp->is_valid_seg(l, r, stp)) {
                auto subamn = getSubset(amns, l, r, stp);
                assert(subamn.size() == part_len);
                double scr;
                if (score->getScoreForSegment(subamn, *prediction, part_id, scr)) {
                    int bg = l;
                    if (stp == -1) {
                        bg = r;
                    }

                    segs.push_back(Segment(bg, bg + part_len - 1, -1, (stp==-1), scr));
                }
            }
        }
    }

    return segs;
}

std::vector<matcher::aacid>
matcher::Matcher::getSubset(std::vector<matcher::aacid> amns, int l, int r, int stp) const {
    std::vector<matcher::aacid> subamn;
    for (int i = l; i != r; i = (i + stp + amns.size())%amns.size()) {
        subamn.push_back(amns[i]);
    }
    subamn.push_back(amns[r]);
    return subamn;
}

int matcher::Matcher::addSegments(const std::vector<matcher::Segment> &part_seg, std::vector<matcher::Segment> &segments,
                                  bool skip_first, bool skip_last, int id) const {

    int cnt_add = 0;
    for (int j = 0; j < part_seg.size(); ++j) {
        if (part_seg[j].l >= skip_first && part_seg[j].r < nrp->getLen() - skip_last) {
            segments.push_back(Segment(part_seg[j].l - skip_first, part_seg[j].r - skip_first, id, part_seg[j].rev, part_seg[j].scor));
            ++cnt_add;
        }
    }

    return cnt_add;
}

void matcher::Matcher::matchSingleUnits(matcher::MatcherBase::Match &match, std::vector<bool> &used_pos) const {
    const double INF = 1e9;
    auto short_parts = prediction->getShortParts();
    int cnt_big_parts = prediction->getNrpsParts().size();
    int n = nrp->getLen();
    int m = short_parts.size();
    bool isSwap = false;
    if (n > m) {
        isSwap = true;
        std::swap(n, m);
    }
    std::vector<std::vector<double>> a(n + 1,
                                       std::vector<double>(m + 1, 1));

    for (int i = 1; i <= n; ++i) {
        for (int j = 1; j <= m; ++j) {
            if (!isSwap) {
                a[i][j] = -score->singleUnitScore(short_parts[j - 1].getAminoacidsPrediction()[0],
                                                 nrp->getAminoacid(i - 1));
                if (used_pos[i - 1] == true) {
                    a[i][j] = 1;
                }
            } else {
                a[i][j] = -score->singleUnitScore(short_parts[i - 1].getAminoacidsPrediction()[0],
                                                 nrp->getAminoacid(j - 1));
                if (used_pos[j - 1] == true) {
                    a[i][j] = 1;
                }
            }
        }
    }

    std::vector<double> u(n + 1), v(m + 1);
    std::vector<int> p(m + 1), way(m + 1);

    for (int i = 1; i <= n; ++i) {
        p[0] = i;
        int j0 = 0;
        std::vector<double> minv(m + 1, INF);
        std::vector<bool> used(m + 1, false);

        do {
            used[j0] = true;
            int i0 = p[j0], delta = INF, j1;
            for (int j=1; j <= m; ++j) {
                if (!used[j]) {
                    int cur = a[i0][j]-u[i0]-v[j];
                    if (cur < minv[j]) {
                        minv[j] = cur, way[j] = j0;
                    }
                    if (minv[j] < delta) {
                        delta = minv[j], j1 = j;
                    }
                }
            }

            for (int j = 0; j <= m; ++j) {
                if (used[j]) {
                    u[p[j]] += delta;
                    v[j] -= delta;
                } else {
                    minv[j] -= delta;
                }
            }
            j0 = j1;
        } while (p[j0] != 0);
        do {
            int j1 = way[j0];
            p[j0] = p[j1];
            j0 = j1;
        } while (j0);
    }

    double newScore = match.score();
    for (int j = 1; j <= m; ++j) {
        int npos = p[j];
        int part_pos = j;
        if (npos == 0) continue;
        if (isSwap) {
            std::swap(npos, part_pos);
        }
        if (used_pos[npos - 1] == false) {
            double curs = score->singleUnitScore(short_parts[part_pos - 1].getAminoacidsPrediction()[0],
                                   nrp->getAminoacid(npos - 1));
            if (curs >= 0) {
                match.match(npos - 1, cnt_big_parts + part_pos - 1, 0);
                newScore += curs;
            }
        }
    }

    match.setScore(newScore);
}

matcher::MatcherBase::Match
matcher::Matcher::getMatch(std::shared_ptr<nrp::NRP> nrp, const nrpsprediction::NRPsPrediction *prediction,
                           const matcher::Score *score) {
    this->nrp = nrp;
    this->prediction = prediction;
    this->score = score;
    return getMatch();
}
