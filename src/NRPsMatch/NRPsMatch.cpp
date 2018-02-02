#include <iostream>
#include <fstream>
#include "NRPsMatch.h"

void nrpsmatch::NRPsMatch::match(int pos, int part_id, int part_pos) {
    parts_id[pos] = part_id;
    parts_pos[pos] = part_pos;
}

int nrpsmatch::NRPsMatch::score() {
    int cnt = 0;
    for (int i = 0; i < parts_id.size(); ++i) {
        if (parts_id[i] >= 0) {
            ++cnt;
        }
    }
    return cnt;
}

std::vector<std::pair<int, int> > nrpsmatch::NRPsMatch::getMatchs() {
    std::vector<std::pair<int, int> > res;
    for (int i = 0; i < parts_id.size(); ++i) {
        res.push_back(std::make_pair(parts_id[i], parts_pos[i]));
    }
    return res;
}

void nrpsmatch::NRPsMatch::print(std::ofstream &out) {
    out << nrp.get_file_name() << "\n";
    int scr = score();
    if (scr < nrp.getLen() - 1) {
        out << "FAIL\n";
    } else {
        out << "MATCH\n";
    }

    out << "number of components : " << nrp.getLen() << "\n";
    for (int i = 0; i < parts_id.size(); ++i) {
        int ri = nrp.getInd(i);
        std::string formula = nrp.getFormula(i);

        out << formula << " -> ";
        if (parts_id[ri] == -1) {
            out << "-\n";
        } else {
            nrpsprediction::AminoacidPrediction amn_pred = nrpParts[parts_id[ri]].getAminoacidsPrediction()[parts_pos[ri]];
            nrpsprediction::AminoacidPrediction::AminoacidProb amprob = amn_pred.getAminoacid(nrp.getAminoacid(ri));
            out << aminoacid::Aminoacids::AMINOACID_NAMES[amprob.aminoacid] << " " << nrpParts[parts_id[ri]].get_orf_name() << " " << parts_pos[ri] << "\n";
        }

    }

    out << nrp.getGraphInString();

    out << "\n\n\n";
}
