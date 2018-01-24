#include <assert.h>
#include <sstream>
#include <iostream>
#include "NRPsPart.h"

std::string nrpsprediction::NRPsPart::get_orf_name() {
    return orf;
}

void nrpsprediction::NRPsPart::add_prediction(int num, std::string predict_aminoacids) {
    assert(num - 1 == aminoacids.size());
    aminoacids.push_back(AminoacidPrediction(num, predict_aminoacids));

}

std::vector<nrpsprediction::AminoacidPrediction> nrpsprediction::NRPsPart::getAminoacidsPrediction() {
    return aminoacids;
}

nrpsprediction::NRPsPart::NRPsPart(std::string file_name, std::string orf_name, int num,
                                   std::string predict_aminoacids) {
    this->file_name = file_name;
    this->orf = orf_name;
    aminoacids.push_back(AminoacidPrediction(num, predict_aminoacids));
}

nrpsprediction::NRPsPart::NRPsPart(std::string file_name, std::string orf_name) {
    this->file_name = file_name;
    this->orf = orf_name;
}
