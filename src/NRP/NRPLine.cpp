#include "NRPLine.h"

nrp::NRPLine::NRPLine(const std::string &file_name, const std::vector<std::string> &strformula,
                      const std::vector<aminoacid::Aminoacids::Aminoacid> &aminoacids, const std::vector<int> &position,
                      const std::string &graph) : NRP(file_name, strformula, aminoacids, position, graph) {}

nrp::NRP::Match nrp::NRPLine::isCover(nrpsprediction::NRPsPrediction nrPsPrediction) {
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


    Match curMatch = isCoverLine(segments, nrPsPrediction, toSmallId, toBigId);


    //TODO max score
    if (curMatch.score() > len - 2) {
        return curMatch;
    }

    Match resMatch(this, nrpparts);
    return resMatch;
}
