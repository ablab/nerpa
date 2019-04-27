//
// Created by olga on 29.01.19.
//

#include <sstream>
#include "PredictionBuilderBase.h"

namespace nrpsprediction {
    const double PredictionBuilderBase::EPS = 1e-4;

    std::pair<std::string, int> PredictionBuilderBase::get_orf_name_and_order(std::string orf) {
        std::stringstream ss(orf);
        std::string prefix;
        std::string orfname;
        std::string tail;
        getline(ss, prefix, '_');
        getline(ss, orfname, '_');
        getline(ss, tail, '_');

        std::stringstream ss2(tail.substr(1));
        int pos;
        ss2 >> pos;

        return std::make_pair(prefix + "_" + orfname, pos);
    }

    int PredictionBuilderBase::getAAbyName(std::string s) {
        for (int i = 0; i < aminoacid::AminoacidInfo::AMINOACID_CNT; ++i) {
            if (s == aminoacid::AminoacidInfo::AMINOACID_NAMES[i]) {
                return i;
            }
        }

        return aminoacid::AminoacidInfo::AMINOACID_CNT - 1;
    }

    NRPsPrediction PredictionBuilderBase::getPrediction() {
        return NRPsPrediction(nrpparts);
    }
}