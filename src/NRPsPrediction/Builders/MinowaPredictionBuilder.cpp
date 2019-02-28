//
// Created by olga on 23.02.19.
//

#include <fstream>
#include <sstream>
#include "MinowaPredictionBuilder.h"

namespace nrpsprediction {
    NRPsPrediction MinowaPredictionBuilder::getPrediction() {
        return NRPsPrediction(nrpparts);
    }

    void MinowaPredictionBuilder::read_file(std::string file_name) {
        std::ifstream in(file_name);
        std::string s;
        while (getline(in, s) && s != "\\") {
            std::string orf_name;
            getline(in, orf_name);
            getline(in, s);
            while(true){

            }
        }

        if (nrpparts.size() > 0 && nrpparts[nrpparts.size() - 1].getAminoacidsPrediction().size() < 2) {
            nrpparts.pop_back();
        }
        in.close();
    }
}
