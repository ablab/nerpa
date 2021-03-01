#include <iostream>
#include <fstream>
#include <algorithm>
#include <Logger/logger.hpp>
#include "NRP/NRP.h"
#include "MatcherBase.h"


void matcher::MatcherBase::Match::match_align(int pos, int part_id, int part_pos) {
    alignment.emplace_back(pos, part_id, part_pos);
}

void matcher::MatcherBase::Match::match(int pos, int part_id, int part_pos) {
    parts_id[pos] = part_id;
    parts_pos[pos] = part_pos;
}

double matcher::MatcherBase::Match::score() const {
    return scr;
}

std::vector<std::pair<int, int> > matcher::MatcherBase::Match::getMatchs() {
    std::vector<std::pair<int, int> > res;
    for (int i = 0; i < nrp->getLen(); ++i) {
        res.push_back(std::make_pair(parts_id[i], parts_pos[i]));
    }
    return res;
}


bool matcher::MatcherBase::Match::isMatched(int i) {
    if (parts_id[i] == -1) {
        return false;
    } else {
        return true;
    }
}


void matcher::MatcherBase::Match::print(std::ostream &out) {
    out << nrp->get_file_name() << " " << nrp->get_extra_info() << "\n";
    if (nrpParts.size() > 0) {
        out << nrpParts[0].get_file_name() << "\n";
    }
    double scr = score();
    out << "SCORE: " << scr << "\n";

    out << "ALIGNMENT:\n";
    out << "ORF_ID A_domain_Idx Prediction_DL-config Prediction_Residue Matched_Residue(Nerpa_score;start_rank-end_rank)";
    out << " -> Monomer_Idx Monomer_Code Monomer_DL-config Monomer_Residue\n";
    for (auto it = alignment.rbegin(); it != alignment.rend(); ++it) {
        int nrp_pos, part_id, part_pos;
        std::tie(nrp_pos, part_id, part_pos) = *it;

        if (part_id == -1) {
            out << "- - - - -";
        } else {
            nrpsprediction::AAdomainPrediction amn_pred = nrpParts[part_id].getAAdomainPrediction()[part_pos];
            nrpsprediction::AAdomainPrediction::AminoacidProb amprob;
            aminoacid::Aminoacid aa = amn_pred.getAAPrediction()[0].aminoacid;
            std::pair<int, int> pos;
            std::pair<double, aminoacid::Aminoacid> res;

            out << nrpParts[part_id].get_orf_name() << " " << part_pos << " ";
            out << aa.getConfiguration() << " " << aa;
            for (auto &mod : amn_pred.get_modifications()) {
                out << "+" << aminoacid::ModificationInfo::NAMES[mod.getId()];
            }
            out << " ";
            if (nrp_pos == -1) {
//                res = std::make_pair(0, amn_pred.getAAPrediction()[0].aminoacid);
                out << "-";
            } else {
                aminoacid::Aminoacid nrp_aa = nrp->getAminoacid(nrp_pos);
                res = scoreFun->getTheBestAAInPred(amn_pred, nrp_aa, amprob, pos);

                out << res.second << "(" << res.first << ";"
                    << pos.first << "-" << pos.second << ")";
            }
        }

        out << " -> ";

        if (nrp_pos == -1) {
            out << "- - - -";
        } else {
            aminoacid::Aminoacid aa = nrp->getAminoacid(nrp_pos);
            out << nrp_pos << " " << nrp->getFormula(nrp_pos) << " " << aa.getConfiguration() << " " << aa.get_name();
            for (auto &mod : aa.getModifications()) {
                out << "+" << aminoacid::ModificationInfo::NAMES[mod.getId()];
            }
        }

        out << "\n";
    }

    out << "GRAPH:\n";
    std::vector<int> rp(parts_id.size());
    for (int i = 0; i < rp.size(); ++i) {
        rp[nrp->getInd(i)] = i;
    }

    out << "number of components : " << parts_id.size() << "\n";
    for (int i = 0; i < parts_id.size(); ++i) {
        int ri = rp[i];
        std::string formula = nrp->getFormula(i);

        out << i << " " << formula << " 0 -> ";

        if (!isMatched(ri)) {
            out << "-\n";
        } else {
            nrpsprediction::AAdomainPrediction amn_pred = nrpParts[parts_id[ri]].getAAdomainPrediction()[parts_pos[ri]];
            nrpsprediction::AAdomainPrediction::AminoacidProb amprob;
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

bool matcher::MatcherBase::Match::operator<(matcher::MatcherBase::Match b) const {
    return this->score() > b.score();
}

void matcher::MatcherBase::Match::print_short(std::ofstream &out) {
    out << nrp->get_file_name() << " " << score() << "; ";
}

void matcher::MatcherBase::Match::print_short_prediction(std::ofstream &out) {
    if (nrpParts.size() == 0) return;

    out << nrpParts[0].get_file_name() << " " << score() << "; ";
}

void matcher::MatcherBase::Match::print_csv(std::ofstream &out) {
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
    int len = parts_id.size();
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

int matcher::MatcherBase::Match::getCntMatch() {
    int cntMatch = 0;
    for (int i = 0; i < parts_id.size(); ++i){
        if (isMatched(i)) {
            cntMatch += 1;
        }
    }
    return cntMatch;
}

void matcher::MatcherBase::Match::setScore(double score) {
    this->scr = score;
}
