#include <iostream>
#include <fstream>
#include <algorithm>
#include "NRP/NRP.h"
#include "Matcher.h"

void matcher::Matcher::Match::match(int pos, int part_id, int part_pos) {
    parts_id[pos] = part_id;
    parts_pos[pos] = part_pos;
}

double matcher::Matcher::Match::score() {
    return scr;
}

std::vector<std::pair<int, int> > matcher::Matcher::Match::getMatchs() {
    std::vector<std::pair<int, int> > res;
    for (int i = 0; i < parts_id.size(); ++i) {
        res.push_back(std::make_pair(parts_id[i], parts_pos[i]));
    }
    return res;
}


bool matcher::Matcher::Match::isMatched(int i) {
    if (parts_id[i] == -1) {
        return false;
    } else {
        return true;
    }
}


void matcher::Matcher::Match::print(std::ofstream &out) {
    out << nrp->get_file_name() << " " << nrp->get_extra_info() << "\n";
    if (nrpParts.size() > 0) {
        out << nrpParts[0].get_file_name() << "\n";
    }
    double scr = score();
    out << "SCORE: " << scr << "("<< nrp->getLen() << ")\n";

    std::vector<int> rp(parts_id.size());
    for (int i = 0; i < rp.size(); ++i) {
        rp[nrp->getInd(i)] = i;
    }

    out << "number of components : " << nrp->getLen() << "\n";
    for (int i = 0; i < parts_id.size(); ++i) {
        int ri = rp[i];
        std::string formula = nrp->getFormula(i); //FIXME i not ri?? check

        out << formula << " -> ";

        if (!isMatched(ri)) {
            out << "-\n";
        } else {
            nrpsprediction::AminoacidPrediction amn_pred = nrpParts[parts_id[ri]].getAminoacidsPrediction()[parts_pos[ri]];


            nrpsprediction::AminoacidPrediction::AminoacidProb amprob;
            std::pair<int, int> pos;
            auto res = scoreFun->getTheBestAAInPred(amn_pred, nrp->getAminoacid(ri), amprob, pos);

            out << res.second << "("  << res.first <<"; "
                << pos.first << "-" << pos.second << ") "
                << nrpParts[parts_id[ri]].get_orf_name() << " " << parts_pos[ri] << "\n";
        }
    }

    out << nrp->getGraphInString();
    out << "\n\n\n";
}

bool matcher::Matcher::Match::operator<(matcher::Matcher::Match b) {
    return this->score() > b.score();
}

void matcher::Matcher::Match::print_short(std::ofstream &out) {
    out << nrp->get_file_name() << " " << score() << "; ";
}

void matcher::Matcher::Match::print_short_prediction(std::ofstream &out) {
    if (nrpParts.size() == 0) return;

    out << nrpParts[0].get_file_name() << " " << score() << "; ";
}

void matcher::Matcher::Match::print_csv(std::ofstream &out) {
    if (nrpParts.size() == 0) {
        return;
    }

    double scr = score();
    out << scr << ",";
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
        if (isMatched(i)) {
            ++cntMatch;
        }
    }

    out << len << ",";
    out << cntMatch << ",";
    out << (len==cntMatch) << ",";
    out << nrp->get_file_name() << ",";
    out << nrpParts[0].get_file_name() << "\n";
}
