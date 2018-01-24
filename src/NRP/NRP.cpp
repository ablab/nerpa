#include <fstream>
#include <iostream>
#include "NRP.h"

void nrp::NRP::parse_fragment_graph(std::string fragment_graph) {
    std::ifstream in(fragment_graph);
    std::string s;
    int line_cnt;
    in >> s >> s >> s >> s >> line_cnt;

    for (int i = 0; i < line_cnt; ++i) {
        int id;
        std::string formula;
        double mass;

        in >> id >> formula >> mass;

        aminoacids.push_back(aminoacid::Aminoacids::get_aminoacid_from_formula(formula));
    }

    in.close();
}

bool nrp::NRP::containNRPsPart(nrpsprediction::NRPsPart predict_part) {
    std::vector<nrpsprediction::AminoacidPrediction> aminoacid_predictions = predict_part.getAminoacidsPrediction();
    for (int i = 0; i < (int)aminoacids.size(); ++i) {
        bool is_ok = true;
        for (int j = 0; j < (int)aminoacid_predictions.size() && is_ok; ++j) {
            int curi = (i + j) % aminoacids.size();
            if (!aminoacid_predictions[j].contain(aminoacids[curi])) {
                is_ok = false;
            }
        }
        if (is_ok == true) {
            return true;
        }
    }

    for (int i = 0; i < (int)aminoacids.size(); ++i) {
        bool is_ok = true;
        int g = 0;
        for (int j = (int)aminoacid_predictions.size() - 1; j >= 0 && is_ok; --j, ++g) {
            int curi = (i + g) % aminoacids.size();
            if (!aminoacid_predictions[j].contain(aminoacids[curi])) {
                is_ok = false;
            }
        }
        if (is_ok == true) {
            return true;
        }
    }


    return false;
}

void nrp::NRP::print() {
    for (int i = 0; i < (int)aminoacids.size(); ++i) {
        std::cerr << aminoacid::Aminoacids::AMINOACID_NAMES[aminoacids[i]] << " ";
    }
    std::cerr << "\n";
}
