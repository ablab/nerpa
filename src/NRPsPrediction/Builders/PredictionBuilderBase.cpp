//
// Created by olga on 29.01.19.
//

#include <sstream>
#include "PredictionBuilderBase.h"

namespace nrpsprediction {
    const double PredictionBuilderBase::EPS = 1e-4;

    std::pair<std::string, int> PredictionBuilderBase::get_orf_name_and_order(std::string orf) {
        std::stringstream ss(orf);

	std::string domain_full_name;
	ss >> domain_full_name;
	int end_of_prefix = domain_full_name.find('_');
	int start_of_tail = domain_full_name.rfind('_');
        std::string prefix = domain_full_name.substr(0, end_of_prefix);
        std::string orfname = domain_full_name.substr(end_of_prefix + 1,
						      start_of_tail - end_of_prefix - 1);
        std::string tail = domain_full_name.substr(start_of_tail);
	
        std::stringstream ss2(tail.substr(1));
        if (tail[0] != 'A' || !(tail[1] >= '0' && tail[1] <= '9')) {
            return std::make_pair(prefix + "_" + orfname, -1);
        }

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

    BgcPrediction PredictionBuilderBase::getPrediction() {
        return BgcPrediction(nrpparts, short_parts);
    }

    PredictionBuilderBase::~PredictionBuilderBase() {

    }
}
