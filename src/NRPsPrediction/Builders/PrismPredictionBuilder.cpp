//
// Created by olga on 01.03.19.
//

#include "PrismPredictionBuilder.h"
#include <Json/json.hpp>
#include <fstream>


nrpsprediction::NRPsPrediction nrpsprediction::PrismPredictionBuilder::getPrediction() {
    return  NRPsPrediction(nrpparts);
}

void nrpsprediction::PrismPredictionBuilder::read_file(std::string file_name) {
    using json = nlohmann::json;

    json jsn;
    std::ifstream in(file_name);
    in >> jsn;

    std::cerr << jsn << "\n";
    in.close();

}
