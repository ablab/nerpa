#include <iostream>
#include "NRPBranch.h"

nrp::NRP::NRPType nrp::NRPBranch::getType() const {
    return NRP::branch_cycle;
}

int nrp::NRPBranch::getLen() const {
    return v1->getLen();
}

int nrp::NRPBranch::getInd(int i) const {
    return v1->getInd(i);
}

std::string nrp::NRPBranch::getFormula(int i) const {
    return v1->getFormula(i);
}

aminoacid::Aminoacid nrp::NRPBranch::getAminoacid(int i) const {
    return v1->getAminoacid(i);
}

void nrp::NRPBranch::print() const {
    v1->print();
}

std::string nrp::NRPBranch::getGraphInString() const {
    return v1->getGraphInString();
}

std::string nrp::NRPBranch::get_file_name() const {
    return v1->get_file_name();
}

std::string nrp::NRPBranch::get_extra_info() const {
    return v1->get_extra_info();
}

std::vector<std::shared_ptr<nrp::NRP>> nrp::NRPBranch::getLines() const {
    std::vector<std::shared_ptr<nrp::NRP>> lines;
    lines.push_back(v1);
    lines.push_back(v2);
    return lines;
}

bool nrp::NRPBranch::is_valid_seg(int l, int r, int stp) const {
    return v1->is_valid_seg(l, r, stp);
}

nrp::NRPBranch::NRPBranch(const nrp::NRP &refNrp) {
    auto lines = refNrp.getLines();
    assert(lines.size() == 2);
    this->v1 = lines[0];
    this->v2 = lines[1];
}

std::string nrp::NRPBranch::structure_to_string() const {
    std::stringstream structure_str;
    int full_len = v1->getFullLen();
    for (int i = (int)full_len - 1; i >= tail_size; --i) {
        structure_str << i;
        if (i != tail_size) {
            structure_str << ",";
        }
    }

    structure_str << "*";
    for (int i = tail_size - 1; i >= 0; --i) {
        structure_str <<"," << i;
    }

    return structure_str.str();
}

int nrp::NRPBranch::getFullLen() const {
    return v1->getFullLen();
}
