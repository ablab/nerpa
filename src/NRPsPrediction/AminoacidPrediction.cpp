#include <sstream>
#include <cmath>
#include <iostream>
#include "AminoacidPrediction.h"


const double nrpsprediction::AminoacidPrediction::EPS = 1e-4;

nrpsprediction::AminoacidPrediction::AminoacidPrediction(int pos, std::string predict_aminoacids) {
    std::stringstream ss(predict_aminoacids);

    std::vector<std::string> tokens;

    while(ss.good()) {
        std::string token;
        getline(ss, token, ';');
        tokens.push_back(token);
    }

    if (tokens[tokens.size() - 1] == "") {
        tokens.pop_back();
    }

    std::vector<std::pair<std::string, double> > aacids;

    for (int i = 0; i < (int)tokens.size(); ++i) {
        aacids.push_back(parse_token(tokens[i]));
    }

    double val = std::max(aacids[2].second, 60.);

    for (int i = 0; i < aacids.size(); ++i) {
        if (aacids[i].second >= val - EPS) {
            aminoacid_prediction.push_back(AminoacidProb(aacids[i].first, aacids[i].second));
            if (aminoacid_prediction[aminoacid_prediction.size() - 1].aminoacid == aminoacid::Aminoacids::Aminoacid::glu) {
                aminoacid_prediction.push_back(AminoacidProb(aminoacid::Aminoacids::Aminoacid::me3_glu,
                                                             aminoacid_prediction[aminoacid_prediction.size() - 1].prob));
            }
            if (aminoacid_prediction[aminoacid_prediction.size() - 1].aminoacid == aminoacid::Aminoacids::Aminoacid::asn) {
                aminoacid_prediction.push_back(AminoacidProb(aminoacid::Aminoacids::Aminoacid::OH_asn,
                                                             aminoacid_prediction[aminoacid_prediction.size() - 1].prob));
            }
        }
    }
}

std::pair<std::string, double> nrpsprediction::AminoacidPrediction::parse_token(std::string token) {
    token.resize(token.size() - 1);
    int pos = token.find('(');
    std::string name = token.substr(0, pos);
    std::string prob = token.substr(pos + 1);

    std::stringstream ss(prob);
    double prb;
    ss >> prb;

    return make_pair(name, prb);
}

bool nrpsprediction::AminoacidPrediction::contain(aminoacid::Aminoacids::Aminoacid aminoacid) {
    for (int i = 0; i < (int)aminoacid_prediction.size(); ++i) {
        if (aminoacid::Aminoacids::same(aminoacid_prediction[i].aminoacid, aminoacid)) {
            return true;
        }
    }
    return false;
}

nrpsprediction::AminoacidPrediction::AminoacidProb
nrpsprediction::AminoacidPrediction::getAminoacid(aminoacid::Aminoacids::Aminoacid aminoacid) const {
    for (int i = 0; i < (int)aminoacid_prediction.size(); ++i) {
        if (aminoacid::Aminoacids::same(aminoacid_prediction[i].aminoacid, aminoacid)) {
            return aminoacid_prediction[i];
        }
    }

    return AminoacidProb(aminoacid::Aminoacids::AMINOACID_CNT, 0);
}

std::pair<int, int> nrpsprediction::AminoacidPrediction::getAmnAcidPos(aminoacid::Aminoacids::Aminoacid aminoacid) const {
    double  prb = -1;
    for (int i = 0; i < (int)aminoacid_prediction.size(); ++i) {
        if (aminoacid::Aminoacids::same(aminoacid_prediction[i].aminoacid, aminoacid)) {
            prb = aminoacid_prediction[i].prob;
            break;
        }
    }
    int bg = -1;
    int ed = -1;

    for (int i = 0; i < (int)aminoacid_prediction.size(); ++i) {
        if (fabs(prb - aminoacid_prediction[i].prob) < EPS) {
            if (bg == -1) {
                bg = i;
            }
            ed = i;
        }
    }

    return std::make_pair(bg, ed);
}
