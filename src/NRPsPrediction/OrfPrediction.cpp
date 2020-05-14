#include <assert.h>
#include <sstream>
#include <iostream>
#include "OrfPrediction.h"

std::string nrpsprediction::OrfPrediction::get_orf_name() const {
    return orf;
}

std::vector<nrpsprediction::AAdomainPrediction> nrpsprediction::OrfPrediction::getAAdomainPrediction() const {
    return aminoacids;
}

nrpsprediction::OrfPrediction::OrfPrediction(std::string file_name, std::string orf_name, bool is_repeatable) {
    this->file_name = file_name;
    this->orf = orf_name;
    this->is_repeatable = is_repeatable;
}

std::string nrpsprediction::OrfPrediction::get_file_name() const {
    return file_name;
}

nrpsprediction::OrfPrediction::OrfPrediction(std::string file_name, std::string orf_name, int num,
                                             const nrpsprediction::AAdomainPrediction &prediction,
                                             bool is_repeatable) {
    this->file_name = file_name;
    this->orf = orf_name;
    this->is_repeatable = is_repeatable;
    aminoacids.push_back(prediction);
}

void nrpsprediction::OrfPrediction::add_prediction(int num, const nrpsprediction::AAdomainPrediction &prediction) {
    assert(num - 1 == aminoacids.size());
    aminoacids.push_back(prediction);
}

nrpsprediction::OrfPrediction::OrfPrediction() {
}
