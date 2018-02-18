#include <iostream>
#include <fstream>
#include <algorithm>
#include "NRP.h"

void nrp::NRP::Match::match(int pos, int part_id, int part_pos) {
    parts_id[pos] = part_id;
    parts_pos[pos] = part_pos;
}

int nrp::NRP::Match::score() {
    int cnt = parts_id.size();
    std::vector<int> difparts;
    for (int i = 0; i < parts_id.size(); ++i) {
        if (parts_id[i] < 0) {
            cnt -= 2;
        } else {
            difparts.push_back(parts_id[i]);
        }
    }

    std::sort(difparts.begin(), difparts.end());
    difparts.resize(std::unique(difparts.begin(), difparts.end()) - difparts.begin());
    cnt -= difparts.size();
    return cnt;
}

std::vector<std::pair<int, int> > nrp::NRP::Match::getMatchs() {
    std::vector<std::pair<int, int> > res;
    for (int i = 0; i < parts_id.size(); ++i) {
        res.push_back(std::make_pair(parts_id[i], parts_pos[i]));
    }
    return res;
}

void nrp::NRP::Match::print(std::ofstream &out) {
    out << nrp->get_file_name() << " " << nrp->get_extra_info() << "\n";
    int scr = score();
    out << "SCORE: " << scr << "("<< nrp->getLen() << ")\n";

    std::vector<int> rp(parts_id.size());
    for (int i = 0; i < rp.size(); ++i) {
        rp[nrp->getInd(i)] = i;
    }

    out << "number of components : " << nrp->getLen() << "\n";
    for (int i = 0; i < parts_id.size(); ++i) {
        int ri = rp[i];
        std::string formula = nrp->getFormula(i);

        out << formula << " -> ";

        if (parts_id[ri] == -1) {
            out << "-\n";
        } else {
            nrpsprediction::AminoacidPrediction amn_pred = nrpParts[parts_id[ri]].getAminoacidsPrediction()[parts_pos[ri]];
            nrpsprediction::AminoacidPrediction::AminoacidProb amprob = amn_pred.getAminoacid(nrp->getAminoacid(ri));
            out << aminoacid::Aminoacids::AMINOACID_NAMES[amprob.aminoacid] << " "
                << nrpParts[parts_id[ri]].get_orf_name() << " " << parts_pos[ri] << "\n";
        }
    }

    out << nrp->getGraphInString();
    out << "\n\n\n";
}

bool nrp::NRP::Match::operator<(nrp::NRP::Match b) {
    return this->score() > b.score();
}

void nrp::NRP::Match::print_short(std::ofstream &out) {
    out << nrp->get_file_name() << " ";
    int scr = score();
    out << scr << "("<< nrp->getLen() << "); ";
}
