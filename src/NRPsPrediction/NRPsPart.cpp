#include <assert.h>
#include <sstream>
#include <iostream>
#include "NRPsPart.h"

std::string nrpsprediction::NRPsPart::get_orf_name() {
    return orf;
}

std::vector<nrpsprediction::AminoacidPrediction> nrpsprediction::NRPsPart::getAminoacidsPrediction() const {
    return aminoacids;
}

nrpsprediction::NRPsPart::NRPsPart(std::string file_name, std::string orf_name) {
    this->file_name = file_name;
    this->orf = orf_name;
}

std::string nrpsprediction::NRPsPart::get_file_name() {
    return file_name;
}

nrpsprediction::NRPsPart::NRPsPart(std::string file_name, std::string orf_name, int num,
                                   const nrpsprediction::AminoacidPrediction &prediction) {
    this->file_name = file_name;
    this->orf = orf_name;
    aminoacids.push_back(prediction);
}

void nrpsprediction::NRPsPart::add_prediction(int num, const nrpsprediction::AminoacidPrediction &prediction) {
    std::cerr << num << " " << aminoacids.size() << "\n";
    assert(num - 1 == aminoacids.size());
    aminoacids.push_back(prediction);
}
