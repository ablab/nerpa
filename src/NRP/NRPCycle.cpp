#include <iostream>
#include "NRPCycle.h"

nrp::NRPCycle::NRPCycle(const std::string &file_name, const std::vector<std::string> &strformula,
                        const std::vector<aminoacid::Aminoacids::Aminoacid> &aminoacids,
                        const std::vector<int> &position, const std::string &graph) : NRP(file_name, strformula,
                                                                                          aminoacids, position,
                                                                                          graph) {
}

nrp::NRP::Match nrp::NRPCycle::isCover(nrpsprediction::NRPsPrediction nrPsPrediction) {
    std::vector<Segment> segments;
    auto nrpparts = nrPsPrediction.getNrpsParts();
    std::vector<int> toSmallId(nrpparts.size(), -1);
    std::vector<int> toBigId;
    int len = this->getLen();

    for (int i = 0; i < nrpparts.size(); ++i) {
        std::vector<nrp::NRP::Segment> part_seg = containNRPsPart(nrpparts[i]);
        for (int j = 0; j < part_seg.size(); ++j) {
            segments.push_back(Segment(part_seg[j].l, part_seg[j].r, i, part_seg[j].rev));
        }

        if (part_seg.size() >= 2) {
            toSmallId[i] = toBigId.size();
            toBigId.push_back(i);
        }
    }

    for (int bg = 0; bg < len; ++bg) {
        std::vector<Segment> tmpSeg;
        for (int i = 0; i < segments.size(); ++i) {
            tmpSeg.push_back(Segment(segments[i].l - bg, segments[i].r - bg, segments[i].part_id, segments[i].rev));
            if (tmpSeg[tmpSeg.size() - 1].l < 0) {
                tmpSeg[tmpSeg.size() - 1].l += len;
                tmpSeg[tmpSeg.size() - 1].r += len;
            }
        }

        Match curMatch = isCoverLine(tmpSeg, nrPsPrediction, toSmallId, toBigId);

        //TODO max score
        if (curMatch.score() > len - 2) {
            std::vector<std::pair<int, int> > matchs = curMatch.getMatchs();
            Match resMatchs(this, nrpparts);
            for (int i = 0; i < matchs.size(); ++i) {
                if (matchs[i].first != -1) {
                    resMatchs.match((i + bg + len) % len, matchs[i].first, matchs[i].second);
                }
            }

            return resMatchs;
        }
    }


    Match resMatch(this, nrpparts);
    return resMatch;
}