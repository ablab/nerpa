#include <assert.h>
#include <sstream>
#include "NRPsPart.h"

std::string nrpsprediction::NRPsPart::get_orf_name() {
    return orf;
}

void nrpsprediction::NRPsPart::add_prediction(int num, std::string predict_aminoacids) {
    assert(num - 1 == aminoacids.size());
    aminoacids.push_back(AminoacidPrediction(num, predict_aminoacids));

}
