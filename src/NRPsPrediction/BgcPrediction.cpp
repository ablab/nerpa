#include <fstream>
#include <sstream>
#include <iostream>
#include <algorithm>
#include <Logger/logger.hpp>
#include "BgcPrediction.h"

namespace nrpsprediction {
    const std::vector<OrfPrediction>& BgcPrediction::getOrfs() const {
        return orfs_;
    }

    const std::vector<OrfPrediction>& BgcPrediction::getShortOrfs() const {
        return short_orfs_;
    }

    int BgcPrediction::getSumPredictionLen() const {
        int len = 0;
        for (int i = 0; i < orfs_.size(); ++i) {
            len += orfs_[i].getAAdomainPrediction().size();
        }

        return len;
    }
}
