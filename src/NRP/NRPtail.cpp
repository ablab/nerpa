#include <iostream>
#include "NRPtail.h"

nrp::NRP::NRPType nrp::NRPtail::getType() const {
    return NRP::branch_cycle;
}

int nrp::NRPtail::getLen() const {
    return v1->getLen();
}

int nrp::NRPtail::getInd(int i) const {
    return v1->getInd(i);
}

std::string nrp::NRPtail::getFormula(int i) const {
    return v1->getFormula(i);
}

aminoacid::Aminoacid nrp::NRPtail::getAminoacid(int i) const {
    return v1->getAminoacid(i);
}

void nrp::NRPtail::print() const {
    v1->print();
}

std::string nrp::NRPtail::getGraphInString() const {
    return v1->getGraphInString();
}

std::string nrp::NRPtail::get_file_name() const {
    return v1->get_file_name();
}

std::string nrp::NRPtail::get_extra_info() const {
    return v1->get_extra_info();
}

std::vector<std::shared_ptr<nrp::NRP>> nrp::NRPtail::getLines() const {
    std::vector<std::shared_ptr<nrp::NRP>> lines;
    lines.push_back(v1);
    lines.push_back(v2);
    return lines;
}

bool nrp::NRPtail::is_valid_seg(int l, int r, int stp) const {
    return v1->is_valid_seg(l, r, stp);
}

nrp::NRPtail::NRPtail(const nrp::NRP &refNrp) {
    auto lines = refNrp.getLines();
    assert(lines.size() == 2);
    this->v1 = lines[0];
    this->v2 = lines[1];
}
