#include <iostream>
#include <fstream>
#include <algorithm>
#include "NRP.h"

void nrp::NRP::Match::match(int pos, int part_id, int part_pos) {
    parts_id[pos] = part_id;
    parts_pos[pos] = part_pos;
}

double nrp::NRP::Match::score() {
    double cnt = 0;
    std::vector<int> difparts;
    for (int i = 0; i < parts_id.size(); ++i) {
        if (parts_id[i] < 0 && (i == 0 || parts_id[i - 1] >= 0)) {
            cnt -= 1;
        } else if (parts_id[i] >= 0) {
            nrpsprediction::AminoacidPrediction amn_pred = nrpParts[parts_id[i]].getAminoacidsPrediction()[parts_pos[i]];
            cnt += amn_pred.getScore(nrp->getAminoacid(i));
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

void nrp::NRP::Match::print(std::ofstream &out, double normScore,  double p_value) {
    out << nrp->get_file_name() << " " << nrp->get_extra_info() << "\n";
    if (nrpParts.size() > 0) {
        out << nrpParts[0].get_file_name() << "\n";
    }
    double scr = score();
    out << "SCORE: " << scr << "("<< nrp->getLen() << ")\n";
    out << "NORMALIZE SCORE: " << normScore << "\n";
    out << "P-VALUE: " << p_value << "\n";

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
            std::pair<int, int> pos = amn_pred.getAmnAcidPos((nrp->getAminoacid(ri)));
            out << aminoacid::Aminoacids::AMINOACID_NAMES[amprob.aminoacid] << "("  << amprob.prob <<"; "
                << pos.first << "-" << pos.second << ") "
                << nrpParts[parts_id[ri]].get_orf_name() << " " << parts_pos[ri] << "\n";
        }
    }

    out << nrp->getGraphInString();
    out << "\n\n\n";
}

bool nrp::NRP::Match::operator<(nrp::NRP::Match b) {
    return this->score() > b.score();
}

void nrp::NRP::Match::print_short(std::ofstream &out, double normScore) {
    out << nrp->get_file_name() << " ";
    out << normScore << "; ";
}

void nrp::NRP::Match::print_short_prediction(std::ofstream &out, double normScore) {
    if (nrpParts.size() == 0) return;

    out << nrpParts[0].get_file_name() << " ";
    out << normScore << "; ";
}

void nrp::NRP::Match::print_csv(std::ofstream &out, double normScore, double p_value) {
    if (nrpParts.size() == 0) {
        return;
    }

    double scr = score();
    out << scr << "," << normScore << ",";
    std::string org = nrp->get_extra_info().substr(0, nrp->get_extra_info().find(' ', 1));
    for (char &j : org) {
        if (j == ',') {
            j = ' ';
        }
    }
    out << org << ",";
    int len = nrp->getLen();
    int cntMatch = 0;
    for (int i = 0; i < parts_id.size(); ++i) {
        std::string formula = nrp->getFormula(i);
        int wasN = 0;
        for (int j = 0; j < formula.size(); ++j) {
            if (formula[j] == 'N') {
                wasN = 1;
            }
        }
        if (wasN == 0) {
            --len;
        }
        if (parts_id[i] != -1) {
            ++cntMatch;
        }
    }

    out << len << ",";
    out << cntMatch << ",";
    out << (len==cntMatch) << ",";
    out << nrp->get_file_name() << ",";
    out << nrpParts[0].get_file_name() << ",";
    out << p_value << "\n";
}