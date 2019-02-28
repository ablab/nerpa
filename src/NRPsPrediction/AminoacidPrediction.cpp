#include <sstream>
#include <cmath>
#include <iostream>
#include <assert.h>
#include "AminoacidPrediction.h"


const double nrpsprediction::AminoacidPrediction::EPS = 1e-4;

bool nrpsprediction::AminoacidPrediction::contain(aminoacid::Aminoacid aminoacid) const {
    for (int i = 0; i < (int)aminoacid_prediction.size(); ++i) {
        if (aminoacid_prediction[i].aminoacid == aminoacid) {
            return true;
        }
    }
    return false;
}

nrpsprediction::AminoacidPrediction::AminoacidProb
nrpsprediction::AminoacidPrediction::getAminoacid(aminoacid::Aminoacid aminoacid) const {
    for (int i = 0; i < (int)aminoacid_prediction.size(); ++i) {
        if (aminoacid_prediction[i].aminoacid == aminoacid) {
            return aminoacid_prediction[i];
        }
    }

    return AminoacidProb(aminoacid::Aminoacid("none"), 0);
}

std::pair<int, int> nrpsprediction::AminoacidPrediction::getAmnAcidPos(aminoacid::Aminoacid aminoacid) const {
    double  prb = -1;
    for (int i = 0; i < (int)aminoacid_prediction.size(); ++i) {
        if (aminoacid_prediction[i].aminoacid == aminoacid) {
            prb = aminoacid_prediction[i].prob;
            break;
        }
    }
    int bg = -1;
    int ed = -1;

    for (int i = 0; i < (int)aminoacid_prediction.size(); ++i) {
        if (fabs(prb - aminoacid_prediction[i].prob) < EPS) {
            if (bg == -1) {
                bg = i;
            }
            ed = i;
        }
    }

    return std::make_pair(bg, ed);
}
