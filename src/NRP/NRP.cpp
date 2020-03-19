#include <fstream>
#include <iostream>
#include <sstream>
#include <algorithm>
#include "NRP.h"

std::string nrp::NRP::getFormula(int i) const {
    return strformula[i];
}

std::string nrp::NRP::getGraphInString() const {
    return graph;
}

int nrp::NRP::getLen() const {
    return len;
}

int nrp::NRP::getInd(int i) const {
    return position[i];
}

aminoacid::Aminoacid nrp::NRP::getAminoacid(int i) const {
    return aminoacids[i];
}

void nrp::NRP::print() const {
    for (int i = 0; i < (int)aminoacids.size(); ++i) {
        std::cerr << aminoacids[i].get_possible_name() << "(" << position[i] << ") ";
    }
    std::cerr << "\n";
}


std::string nrp::NRP::get_file_name() const {
    return file_name;
}

std::string nrp::NRP::get_extra_info() const {
    return extra_info;
}

std::vector<aminoacid::Aminoacid> nrp::NRP::getAminoacids() const {
    return aminoacids;
}

void nrp::NRP::deleteAA(int i) {
    for (int j = i; j < len - 1; ++j) {
        std::swap(aminoacids[j], aminoacids[j + 1]);
        std::swap(position[j], position[j + 1]);
        //std::swap(strformula[j], strformula[j + 1]);
    }
    len -= 1;
}

int nrp::NRP::getFullLen() const {
    return aminoacids.size();
}

void nrp::NRP::insertAA(int i) {
    aminoacids.push_back(aminoacid::Aminoacid("none"));
    position.push_back(len);
    strformula.push_back("-");

    for (int j = len; j > i; --j) {
        std::swap(aminoacids[j], aminoacids[j - 1]);
        std::swap(position[j], position[j - 1]);
        //std::swap(strformula[j], strformula[j - 1]);
    }

    len += 1;
}

std::string nrp::NRP::structure_to_string() const {
    std::stringstream structure_str;
    for (int i = 0; i < aminoacids.size(); ++i) {
        structure_str << i;
        if (i != (int)aminoacids.size() - 1) {
            structure_str << ",";
        }
    }
    return structure_str.str();
}
