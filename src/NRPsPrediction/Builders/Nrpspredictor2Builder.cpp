//
// Created by olga on 29.01.19.
//

#include <fstream>
#include <sstream>
#include "Nrpspredictor2Builder.h"

namespace nrpsprediction {
    void Nrpspredictor2Builder::read_file(std::string file_name) {
        std::ifstream in(file_name);
        std::string s;
        while (getline(in, s)) {
            std::stringstream ss(s);
            std::string orf_name;
            std::string predict_aminoacid;
            std::string predict_aminoacids;

            ss >> orf_name >> predict_aminoacid >> predict_aminoacids;
            std::pair <std::string, int> orf_name_num = get_orf_name_and_order(orf_name);

            if (!nrpparts.empty() && nrpparts[nrpparts.size() - 1].get_orf_name() == orf_name_num.first) {
                nrpparts[nrpparts.size() - 1].add_prediction(orf_name_num.second,
                                                             AAdomainPrediction(orf_name_num.second, parse_predictions(predict_aminoacids)));
            } else {
                if (nrpparts.size() > 0 && nrpparts[nrpparts.size() - 1].getAAdomainPrediction().size() < 2) {
                    nrpparts.pop_back();
                }
                nrpparts.push_back(OrfPrediction(file_name, orf_name_num.first, orf_name_num.second,
                                                 AAdomainPrediction(orf_name_num.second, parse_predictions(predict_aminoacids))));
            }
        }
        if (nrpparts.size() > 0 && nrpparts[nrpparts.size() - 1].getAAdomainPrediction().size() < 2) {
            nrpparts.pop_back();
        }
        in.close();
    }

    std::vector<AAdomainPrediction::AminoacidProb> Nrpspredictor2Builder::parse_predictions(std::string predictions) {
        std::vector<AAdomainPrediction::AminoacidProb> aminoacid_prediction;
        std::stringstream ss(predictions);

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
                aminoacid_prediction.push_back(
                        AAdomainPrediction::AminoacidProb(
                                aminoacid::Aminoacid(getAAbyName(aacids[i].first)), aacids[i].second));
            }
        }

        return aminoacid_prediction;
    }

    std::pair<std::string, double> Nrpspredictor2Builder::parse_token(std::string token) {
        token.resize(token.size() - 1);
        int pos = token.find('(');
        std::string name = token.substr(0, pos);
        std::string prob = token.substr(pos + 1);

        std::stringstream ss(prob);
        double prb;
        ss >> prb;

        return make_pair(name, prb);
    }
}