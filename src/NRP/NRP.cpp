#include <fstream>
#include <iostream>
#include <sstream>
#include "NRP.h"

void nrp::NRP::parse_fragment_graph(std::string fragment_graph) {
    file_name = fragment_graph;
    std::ifstream in(fragment_graph);
    std::string s;
    int line_cnt;
    in >> s >> s >> s >> s >> line_cnt;

    for (int i = 0; i < line_cnt; ++i) {
        int id;
        std::string formula;
        double mass;

        in >> id >> formula >> mass;

        std::stringstream ss;
        ss << id << " " << formula << " " << mass;

        strformula.push_back(ss.str());

        aminoacids.push_back(aminoacid::Aminoacids::get_aminoacid_from_formula(formula));
    }

    std::string tmp;
    while(getline(in, tmp)) {
        graph += tmp;
        graph += "\n";
    }

    in.close();
}


std::string nrp::NRP::getFormula(int i) {
    return strformula[i];
}

std::string nrp::NRP::getGraphInString() {
    return graph;
}

std::vector<nrp::NRP::Segment> nrp::NRP::containNRPsPart(nrpsprediction::NRPsPart predict_part) {
    std::vector<Segment> res;
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
            res.push_back(Segment(i, i + aminoacid_predictions.size() - 1, 0));
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
            res.push_back(Segment(i, i + aminoacid_predictions.size() - 1, 1));
        }
    }

    return res;
}


int nrp::NRP::getLen() {
    return aminoacids.size();
}

int nrp::NRP::getInd(int i) {
    return i;
}

aminoacid::Aminoacids::Aminoacid nrp::NRP::getAminoacid(int i) {
    return aminoacids[i];
}

void nrp::NRP::print() {
    for (int i = 0; i < (int)aminoacids.size(); ++i) {
        std::cerr << aminoacid::Aminoacids::AMINOACID_NAMES[aminoacids[i]] << " ";
    }
    std::cerr << "\n";
}


std::string nrp::NRP::get_file_name() {
    return file_name;
}
