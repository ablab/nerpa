#include <sstream>
#include "AminoacidPrediction.h"

nrpsprediction::AminoacidPrediction::AminoacidPrediction(int pos, std::string predict_aminoacids) {
    std::stringstream ss(predict_aminoacids);

    std::vector<std::string> tokens;

    while(ss.good()) {
        std::string token;
        getline(ss, token, ';');
        tokens.push_back(token);
    }

    std::vector<std::pair<std::string, double> > aacids;

    for (int i = 0; i < (int)tokens.size(); ++i) {
        aacids.push_back(parse_token(tokens[i]));
    }

    double val = aacids[2].second;

    for (int i = 0; i < aacids.size(); ++i) {
        if (aacids[i].second >= val - EPS) {
            aminoacid_prediction.push_back(AminoacidProb(aacids[i].first, aacids[i].second));
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
