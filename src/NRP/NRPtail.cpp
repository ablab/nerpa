#include <iostream>
#include "NRPtail.h"

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

std::vector<nrp::NRPLine> nrp::NRPtail::getLines() const {
    std::vector<nrp::NRPLine> lines;
    lines.push_back(v1);
    lines.push_back(v2);
    return lines;
}
