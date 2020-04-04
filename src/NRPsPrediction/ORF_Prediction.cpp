#include <assert.h>
#include <sstream>
#include <iostream>
#include "ORF_Prediction.h"

std::string nrpsprediction::ORF_Prediction::get_orf_name() const {
    return orf;
}

std::vector<nrpsprediction::AAdomain_Prediction> nrpsprediction::ORF_Prediction::getAminoacidsPrediction() const {
    return aminoacids;
}

nrpsprediction::ORF_Prediction::ORF_Prediction(std::string file_name, std::string orf_name) {
    this->file_name = file_name;
    this->orf = orf_name;
}

std::string nrpsprediction::ORF_Prediction::get_file_name() const {
    return file_name;
}

nrpsprediction::ORF_Prediction::ORF_Prediction(std::string file_name, std::string orf_name, int num,
                                               const nrpsprediction::AAdomain_Prediction &prediction) {
    this->file_name = file_name;
    this->orf = orf_name;
    aminoacids.push_back(prediction);
}

void nrpsprediction::ORF_Prediction::add_prediction(int num, const nrpsprediction::AAdomain_Prediction &prediction) {
    assert(num - 1 == aminoacids.size());
    aminoacids.push_back(prediction);
}

nrpsprediction::ORF_Prediction::ORF_Prediction() {
}
