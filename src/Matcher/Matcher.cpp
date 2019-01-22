//
// Created by olga on 22.01.19.
//

#include <NRP/NRPLine.h>
#include "Matcher.h"

nrp::NRP::Match matcher::Matcher::getMatch() const {
    if (nrp.getType() == nrp::NRP::line) {
        return getLineMatch();
    } else if (nrp.getType() == nrp::NRP::cycle) {
        return getCycleMatch();
    } else {
        return getBranchMatch();
    }
}

nrp::NRP::Match matcher::Matcher::getLineMatch() const {
    std::vector<nrp::NRP::Segment> segments;
    auto nrpparts = prediction.getNrpsParts();
    std::vector<int> toSmallId(nrpparts.size(), -1);
    std::vector<int> toBigId;
    int len = nrp.getLen();

    for (int i = 0; i < nrpparts.size(); ++i) {
        std::vector<nrp::NRP::Segment> part_seg = nrp.containNRPsPart(nrpparts[i]);
        for (int j = 0; j < part_seg.size(); ++j) {
            segments.push_back(nrp::NRP::Segment(part_seg[j].l, part_seg[j].r, i, part_seg[j].rev, part_seg[j].scor));
        }

        if (part_seg.size() >= 2) {
            toSmallId[i] = toBigId.size();
            toBigId.push_back(i);
        }
    }

    double best_score = -len-1;
    nrp::NRP::Match resMatch(&nrp, nrpparts);

    nrp::NRP::Match curMatch = nrp.isCoverLine(segments, prediction, toSmallId, toBigId);

    if (curMatch.score() > best_score) {
        std::vector<std::pair<int, int> > matchs = curMatch.getMatchs();

        for (int i = 0; i < matchs.size(); ++i) {
            resMatch.match(i, matchs[i].first, matchs[i].second);
        }
    }

    return resMatch;
}

nrp::NRP::Match matcher::Matcher::getCycleMatch() const {
    std::vector<nrp::NRP::Segment> segments;
    auto nrpparts = prediction.getNrpsParts();
    std::vector<int> toSmallId(nrpparts.size(), -1);
    std::vector<int> toBigId;
    int len = nrp.getLen();

    for (int i = 0; i < nrpparts.size(); ++i) {
        std::vector<nrp::NRP::Segment> part_seg = nrp.containNRPsPart(nrpparts[i]);
        for (int j = 0; j < part_seg.size(); ++j) {
            segments.push_back(nrp::NRP::Segment(part_seg[j].l, part_seg[j].r, i, part_seg[j].rev, part_seg[j].scor));
        }

        if (part_seg.size() >= 2) {
            toSmallId[i] = toBigId.size();
            toBigId.push_back(i);
        }
    }

    double best_score = -len - 1;
    nrp::NRP::Match resMatchs(&nrp, nrpparts);

    for (int bg = 0; bg < len; ++bg) {
        std::vector<nrp::NRP::Segment> tmpSeg;
        for (int i = 0; i < segments.size(); ++i) {
            tmpSeg.push_back(nrp::NRP::Segment(segments[i].l - bg, segments[i].r - bg, segments[i].part_id, segments[i].rev, segments[i].scor));
            if (tmpSeg[tmpSeg.size() - 1].l < 0) {
                tmpSeg[tmpSeg.size() - 1].l += len;
                tmpSeg[tmpSeg.size() - 1].r += len;
            }
        }

        nrp::NRP::Match curMatch = updateMatch(prediction, nrp.isCoverLine(tmpSeg, prediction, toSmallId, toBigId), bg);
        if (curMatch.score() > best_score) {
            resMatchs = curMatch;
            best_score = curMatch.score();
        }
    }


    return resMatchs;
}


nrp::NRP::Match matcher::Matcher::getBranchMatch() const {
    std::vector<nrp::NRPLine> lines = nrp.getLines();
    nrp::NRPLine v1 = lines[0];
    nrp::NRPLine v2 = lines[1];

    Matcher matcher1(v1, prediction);
    Matcher matcher2(v2, prediction);
    nrp::NRP::Match m1 = matcher1.getMatch();
    nrp::NRP::Match m2 = matcher2.getMatch();
    if (m1.score() > m2.score()) {
        return m1;
    } else {
        return m2;
    }
}

nrp::NRP::Match
matcher::Matcher::updateMatch(const nrpsprediction::NRPsPrediction &nrPsPrediction, nrp::NRP::Match match,
                              int bg) const {
    nrp::NRP::Match nmatch(&nrp, nrPsPrediction.getNrpsParts());
    std::vector<std::pair<int, int> > part_id_pos = match.getMatchs();
    for (int i = 0; i < part_id_pos.size(); ++i) {
        nmatch.match((i + bg)%nrp.getLen(), part_id_pos[i].first, part_id_pos[i].second);
    }
    return nmatch;
}

