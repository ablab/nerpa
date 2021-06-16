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
    out
            << "ORF_ID A_domain_Idx Prediction_DL-config Prediction_Top_Residue Prediction_Top_Score Prediction_Modifications";
    out << " Matched_Residue Matched_Residue_Score Nerpa_Score";
    out << " Monomer_Idx Monomer_Code Monomer_DL-config Monomer_Residue Monomer_Modifications\n";

    auto _print_modifications = [&out](const std::vector<aminoacid::Modification> &mods) {
        if (mods.empty()) {
            out << "-";
        }
        for (auto &mod : mods) {
            if (&mod != &mods[0]) {
                out << "+";
            }
            out << aminoacid::ModificationInfo::NAMES[mod.getId()];
        }
    };

    for (auto it = alignment.rbegin(); it != alignment.rend(); ++it) {
        int nrp_pos, part_id, part_pos;
        std::tie(nrp_pos, part_id, part_pos) = *it;

        if (part_id == -1) {
            out << "- - - - - - - - " << scoreFun->InsertionScore();
        } else {
            nrpsprediction::AAdomainPrediction amn_pred = nrpParts[part_id].getAAdomainPrediction()[part_pos];
            nrpsprediction::AAdomainPrediction::AminoacidProb amprob;
            nrpsprediction::AAdomainPrediction::AminoacidProb top_amprob = amn_pred.getAAPrediction()[0];
            std::pair<int, int> pos;
            std::pair<double, aminoacid::Aminoacid> res;
            double top_score = top_amprob.prob;

            // ORF_ID A_domain_Idx
            out << nrpParts[part_id].get_orf_name() << " " << part_pos << " ";
            // Prediction_DL-config Prediction_Top_Residue Prediction_Top_Score
            out << top_amprob.aminoacid.getConfiguration() << " " << top_amprob.aminoacid << " " << top_score
                << " ";
            // Prediction_Modifications
            _print_modifications(amn_pred.get_modifications());
            out << " ";
            if (nrp_pos == -1) {
                out << "none none " << scoreFun->DeletionScore();
            } else {
                aminoacid::Aminoacid nrp_aa = nrp->getAminoacid(nrp_pos);
                res = scoreFun->getTheBestAAInPred(amn_pred, nrp_aa, amprob, pos);
                // Matched_Residue Matched_Residue_Score Nerpa_Score
                out << res.second << " " << amprob.prob << " " << scoreFun->aaScore(amn_pred, nrp_aa);
            }
        }

        out << " ";

        if (nrp_pos == -1) {
            out << "- - - - -";
        } else {
            aminoacid::Aminoacid aa = nrp->getAminoacid(nrp_pos);
            // Monomer_Idx Monomer_Code Monomer_DL-config Monomer_Residue
            out << nrp_pos << " " << nrp->getFormula(nrp_pos) << " " << aa.getConfiguration() << " " << aa.get_name()
                << " ";
            // Monomer_Modifications
            _print_modifications(aa.getModifications());
        }

        out << "\n";
    }
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

matcher::MatcherBase::~MatcherBase() {

}
