//
// Created by olga on 29.01.19.
//

#include <sstream>
#include "PredictionBuilderBase.h"

namespace nrpsprediction {
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

        return std::make_pair(orfname, pos);
    }
}