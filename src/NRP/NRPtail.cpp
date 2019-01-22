#include <iostream>
#include "NRPtail.h"

nrp::NRP::Match nrp::NRPtail::isCover(const nrpsprediction::NRPsPrediction& nrPsPrediction) const {
    Match m1 = v1.isCover(nrPsPrediction);
    Match m2 = v2.isCover(nrPsPrediction);
    if (m1.score() > m2.score()) {
        return m1;
    } else {
        return m2;
    }
}

std::vector<nrp::NRP::Segment> nrp::NRPtail::containNRPsPart(nrpsprediction::NRPsPart predict_part) const {
    std::vector<nrp::NRP::Segment> seg1 = v1.containNRPsPart(predict_part), seg2 = v2.containNRPsPart(predict_part);

    seg1.insert(seg1.end(), seg2.begin(), seg2.end());
    return seg1;
}

nrp::NRP::NRPType nrp::NRPtail::getType() const {
    return NRP::branch_cycle;
}

int nrp::NRPtail::getLen() const {
    return v1.getLen();
}

int nrp::NRPtail::getInd(int i) const {
    return v1.getInd(i);
}

std::string nrp::NRPtail::getFormula(int i) const {
    return v1.getFormula(i);
}

aminoacid::Aminoacids::Aminoacid nrp::NRPtail::getAminoacid(int i) const {
    return v1.getAminoacid(i);
}

void nrp::NRPtail::print() const {
    v1.print();
}

std::string nrp::NRPtail::getGraphInString() const {
    return v1.getGraphInString();
}

std::string nrp::NRPtail::get_file_name() const {
    return v1.get_file_name();
}

std::string nrp::NRPtail::get_extra_info() const {
    return v1.get_extra_info();
}
