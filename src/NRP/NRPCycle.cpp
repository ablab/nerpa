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

    int best_score = 0;
    Match resMatchs(this, nrpparts);

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

        if (curMatch.score() > best_score) {
            std::vector<std::pair<int, int> > matchs = curMatch.getMatchs();

            for (int i = 0; i < matchs.size(); ++i) {
                resMatchs.match((i + bg + len) % len, matchs[i].first, matchs[i].second);
            }
        }
    }


    return resMatchs;
}

std::vector<nrp::NRP::Segment> nrp::NRPCycle::containNRPsPart(nrpsprediction::NRPsPart predict_part) {
    std::vector<Segment> res;
    std::vector<nrpsprediction::AminoacidPrediction> aminoacid_predictions = predict_part.getAminoacidsPrediction();
    for (int i = 0; i < (int)aminoacids.size(); ++i) {
        bool is_ok = true;
        for (int j = 0; j < (int)aminoacid_predictions.size() && is_ok; ++j) {
            int curi = (i + j) % aminoacids.size();
            if (!aminoacid_predictions[j].contain(aminoacids[curi])) {
                is_ok = false;
            }
        }
        if (is_ok == true) {
            res.push_back(Segment(i, i + aminoacid_predictions.size() - 1, -1, 0));
        }
    }

    for (int i = 0; i < (int)aminoacids.size(); ++i) {
        bool is_ok = true;
        int g = 0;
        for (int j = (int)aminoacid_predictions.size() - 1; j >= 0 && is_ok; --j, ++g) {
            int curi = (i + g) % aminoacids.size();
            if (!aminoacid_predictions[j].contain(aminoacids[curi])) {
                is_ok = false;
            }
        }
        if (is_ok == true) {
            res.push_back(Segment(i, i + aminoacid_predictions.size() - 1, -1, 1));
        }
    }

    return res;
}
