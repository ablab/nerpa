//
// Created by olga on 22.01.19.
//

#include <NRP/NRPLine.h>
#include <algorithm>
#include "Matcher.h"

matcher::Matcher::Match matcher::Matcher::getMatch() const {
    if (nrp.getType() == nrp::NRP::line) {
        return getLineMatch(true, true);
    } else if (nrp.getType() == nrp::NRP::cycle) {
        return getCycleMatch();
    } else {
        return getBranchMatch();
    }
}

matcher::Matcher::Match matcher::Matcher::getLineMatch(bool can_skip_first, bool can_skip_last) const {
    std::vector<Segment> segments;
    auto nrpparts = prediction.getNrpsParts();
    std::vector<int> toSmallId(nrpparts.size(), -1);
    std::vector<int> toBigId;
    int len = nrp.getLen();
    bool skip_first = (can_skip_first && len > 0 && (nrp.getAminoacid(0).get_name() == "none"));
    bool skip_last = (can_skip_last && len > 0 && (nrp.getAminoacid(len - 1).get_name() == "none"));

    len -= (skip_first + skip_last);

    for (int i = 0; i < nrpparts.size(); ++i) {
        std::vector<Segment> part_seg = matche_seg(nrpparts[i]);

        int cnt_add = addSegments(part_seg, segments, skip_first, skip_last, i);

        if (cnt_add >= 2) {
            toSmallId[i] = toBigId.size();
            toBigId.push_back(i);
        }
    }

    return updateMatch(prediction, isCoverLine(segments, toSmallId, toBigId, len), skip_first);
}

matcher::Matcher::Match matcher::Matcher::getCycleMatch() const {
    std::vector<Segment> segments;
    auto nrpparts = prediction.getNrpsParts();
    std::vector<int> toSmallId(nrpparts.size(), -1);
    std::vector<int> toBigId;
    int len = nrp.getLen();

    for (int i = 0; i < nrpparts.size(); ++i) {
        std::vector<Segment> part_seg = matche_seg(nrpparts[i]);

        for (int j = 0; j < part_seg.size(); ++j) {
            segments.push_back(Segment(part_seg[j].l, part_seg[j].r, i, part_seg[j].rev, part_seg[j].scor));
        }

        if (part_seg.size() >= 2) {
            toSmallId[i] = toBigId.size();
            toBigId.push_back(i);
        }
    }

    double best_score = score.minScore(len);
    matcher::Matcher::Match resMatchs(&nrp, nrpparts, score.minScore(len));

    for (int bg = 0; bg < len; ++bg) {
        std::vector<Segment> tmpSeg;
        for (int i = 0; i < segments.size(); ++i) {
            tmpSeg.push_back(Segment(segments[i].l - bg, segments[i].r - bg, segments[i].part_id, segments[i].rev, segments[i].scor));
            if (tmpSeg[tmpSeg.size() - 1].l < 0) {
                tmpSeg[tmpSeg.size() - 1].l += len;
                tmpSeg[tmpSeg.size() - 1].r += len;
            }
        }

        matcher::Matcher::Match curMatch = updateMatch(prediction, isCoverLine(tmpSeg, toSmallId, toBigId, len), bg);
        if (curMatch.score() > best_score) {
            resMatchs = curMatch;
            best_score = curMatch.score();
        }
    }


    return resMatchs;
}


matcher::Matcher::Match matcher::Matcher::getBranchMatch() const {
    std::vector<nrp::NRPLine> lines = nrp.getLines();
    nrp::NRPLine v1 = lines[0];
    nrp::NRPLine v2 = lines[1];

    assert(v1.getFormula(0) == v2.getFormula(0));

    Matcher matcher1(v1, prediction);
    Matcher matcher2(v2, prediction);
    matcher::Matcher::Match m1 = matcher1.getLineMatch(true, false);
    matcher::Matcher::Match m2 = matcher2.getLineMatch(true, false);
    if (m1.score() > m2.score()) {
        return m1;
    } else {
        return m2;
    }
}

matcher::Matcher::Match
matcher::Matcher::updateMatch(const nrpsprediction::NRPsPrediction &nrPsPrediction, matcher::Matcher::Match match,
                              int bg) const {
    matcher::Matcher::Match nmatch(&nrp, nrPsPrediction.getNrpsParts(), match.score());
    std::vector<std::pair<int, int> > part_id_pos = match.getMatchs();
    for (int i = 0; i < part_id_pos.size(); ++i) {
        nmatch.match((i + bg)%nrp.getLen(), part_id_pos[i].first, part_id_pos[i].second);
    }
    return nmatch;
}

matcher::Matcher::Match
matcher::Matcher::isCoverLine(std::vector<Segment> &segments,
                      const std::vector<int> &toSmallId, const std::vector<int> &toBigId, int len) const {
    std::sort(segments.begin(), segments.end());

    std::vector<std::vector<std::vector<double> > > d(len + 1,
                                                      std::vector<std::vector<double> >((1 << toBigId.size()),
                                                                                        std::vector<double>(2, score.minScore(len))));
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
            if (pos != 0 && d[pos][msk][0] < d[pos - 1][msk][1] + score.openGap()) {
                d[pos][msk][0] = d[pos - 1][msk][1] + score.openGap();
                pgp[pos][msk][0] = 1;
            }
            if (pos != 0 && d[pos][msk][0] < d[pos - 1][msk][0] + score.continueGap()) {
                d[pos][msk][0] = d[pos - 1][msk][0] + score.continueGap();
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
                    if (d[segments[curseg].r + 1][nmsk][1] < d[pos][msk][gp] + score.addSegment(segments[curseg])) {
                        d[segments[curseg].r + 1][nmsk][1] = d[pos][msk][gp] + score.addSegment(segments[curseg]);
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

    double mn = score.minScore(len);
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

    matcher::Matcher::Match nrPsMatch(&nrp, prediction.getNrpsParts(), mn);
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


std::vector<matcher::Segment> matcher::Matcher::matche_seg(const nrpsprediction::NRPsPart &predict_part) const {
    std::vector<Segment> segs;
    std::vector<aacid> amns = nrp.getAminoacids();
    int part_len = predict_part.getAminoacidsPrediction().size();

    if (part_len > amns.size() || part_len == 0) {
        return segs;
    }

    for (int l = 0; l < amns.size(); ++l) {
        for (int stp = -1; stp < 2; stp += 2) {
            int r =  (l + stp*(part_len - 1) + amns.size()) % amns.size();
            if (nrp.is_valid_seg(l, r, stp)) {
                auto subamn = getSubset(amns, l, r, stp);
                assert(subamn.size() == part_len);
                double scr;
                if (score.getScoreForSegment(subamn, predict_part, scr)) {
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
        if (part_seg[j].l >= skip_first && part_seg[j].r < nrp.getLen() - skip_last) {
            segments.push_back(Segment(part_seg[j].l - skip_first, part_seg[j].r - skip_first, id, part_seg[j].rev, part_seg[j].scor));
            ++cnt_add;
        }
    }

    return cnt_add;
}


