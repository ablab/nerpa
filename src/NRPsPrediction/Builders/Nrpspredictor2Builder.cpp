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
                nrpparts[nrpparts.size() - 1].add_prediction(orf_name_num.second, predict_aminoacids);
            } else {
                if (nrpparts.size() > 0 && nrpparts[nrpparts.size() - 1].getAminoacidsPrediction().size() < 2) {
                    nrpparts.pop_back();
                }
                nrpparts.push_back(NRPsPart(file_name, orf_name_num.first, orf_name_num.second, predict_aminoacids));
            }
        }
        if (nrpparts.size() > 0 && nrpparts[nrpparts.size() - 1].getAminoacidsPrediction().size() < 2) {
            nrpparts.pop_back();
        }
        in.close();
    }

    NRPsPrediction Nrpspredictor2Builder::getPrediction() {
        return  NRPsPrediction(nrpparts);
    }
}