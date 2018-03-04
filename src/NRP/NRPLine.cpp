#include "NRPLine.h"

nrp::NRPLine::NRPLine(const std::string &file_name, const std::vector<std::string> &strformula,
                      const std::vector<aminoacid::Aminoacids::Aminoacid> &aminoacids, const std::vector<int> &position,
                      const std::string &graph, const std::string& extra_info) : NRP(file_name, strformula, aminoacids, position, graph, extra_info) {}

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

    int best_score = -len-1;
    Match resMatch(this, nrpparts);

    Match curMatch = isCoverLine(segments, nrPsPrediction, toSmallId, toBigId);

    if (curMatch.score() > best_score) {
        std::vector<std::pair<int, int> > matchs = curMatch.getMatchs();

        for (int i = 0; i < matchs.size(); ++i) {
            resMatch.match(i, matchs[i].first, matchs[i].second);
        }
    }

    return resMatch;
}

std::vector<nrp::NRP::Segment> nrp::NRPLine::containNRPsPart(nrpsprediction::NRPsPart predict_part) {
    std::vector<Segment> res;
    std::vector<nrpsprediction::AminoacidPrediction> aminoacid_predictions = predict_part.getAminoacidsPrediction();
    for (int i = 0; i < (int)aminoacids.size() - (int)aminoacid_predictions.size() + 1; ++i) {
        bool is_ok = true;
        for (int j = 0; j < (int)aminoacid_predictions.size() && is_ok; ++j) {
            int curi = i + j;
            if (!aminoacid_predictions[j].contain(aminoacids[curi])) {
                is_ok = false;
            }
        }
        if (is_ok == true) {
            res.push_back(Segment(i, i + aminoacid_predictions.size() - 1, -1, 0));
        }
    }

    for (int i = 0; i < (int)aminoacids.size() - (int)aminoacid_predictions.size() + 1; ++i) {
        bool is_ok = true;
        int g = 0;
        for (int j = (int)aminoacid_predictions.size() - 1; j >= 0 && is_ok; --j, ++g) {
            int curi = i + g;
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

nrp::NRP::NRPType nrp::NRPLine::getType() {
    return NRP::line;
}
