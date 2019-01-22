#include <iostream>
#include <Logger/logger.hpp>
#include "NRPCycle.h"

nrp::NRPCycle::NRPCycle(const std::string &file_name, const std::vector<std::string> &strformula,
                        const std::vector<aminoacid::Aminoacids::Aminoacid> &aminoacids,
                        const std::vector<int> &position, const std::string &graph, const std::string& extra_info) : NRP(file_name, strformula,
                                                                                          aminoacids, position,
                                                                                          graph, extra_info) {
}

std::vector<nrp::NRP::Segment> nrp::NRPCycle::containNRPsPart(nrpsprediction::NRPsPart predict_part) const {
    std::vector<Segment> res;
    std::vector<nrpsprediction::AminoacidPrediction> aminoacid_predictions = predict_part.getAminoacidsPrediction();
    for (int i = 0; i < (int)aminoacids.size(); ++i) {
        int cnt_mismatch = 0;
        double segscor = 0;
        for (int j = 0; j < (int)aminoacid_predictions.size() && cnt_mismatch < 2; ++j) {
            int curi = (i + j) % aminoacids.size();
            if (!aminoacid_predictions[j].contain(aminoacids[curi])) {
                cnt_mismatch += 1;
            }
            segscor += aminoacid_predictions[j].getScore(aminoacids[curi]);
        }

        if (cnt_mismatch == 0 || (cnt_mismatch == 1 && aminoacid_predictions.size() > 4)) {
            res.push_back(Segment(i, i + aminoacid_predictions.size() - 1, -1, 0, segscor));
        }
    }

    for (int i = 0; i < (int)aminoacids.size(); ++i) {
        int cnt_mismatch = 0;
        int g = 0;
        double segscor = 0;
        for (int j = (int)aminoacid_predictions.size() - 1; j >= 0 && cnt_mismatch < 2; --j, ++g) {
            int curi = (i + g) % aminoacids.size();
            if (!aminoacid_predictions[j].contain(aminoacids[curi])) {
                cnt_mismatch += 1;
            }
            segscor += aminoacid_predictions[j].getScore(aminoacids[curi]);
        }

        if (cnt_mismatch == 0 || (cnt_mismatch == 1 && aminoacid_predictions.size() > 4)) {
            res.push_back(Segment(i, i + aminoacid_predictions.size() - 1, -1, 1, segscor));
        }
    }

    return res;
}

nrp::NRP::NRPType nrp::NRPCycle::getType() const {
    return NRP::cycle;
}

//TODO
std::vector<nrp::NRPLine> nrp::NRPCycle::getLines() const {
    ERROR("Get lines for NRP cycle. NOT implemented");
    return std::vector<nrp::NRPLine>();
}
