#include <fstream>
#include <sstream>
#include <iostream>
#include <algorithm>
#include <Logger/logger.hpp>
#include "NRPsPrediction.h"

namespace nrpsprediction {
    const std::vector<NRPsPart>& NRPsPrediction::getNrpsParts() const {
        return nrpparts;
    }

    const std::vector<NRPsPart>& NRPsPrediction::getShortParts() const {
        return short_parts;
    }

    int NRPsPrediction::getSumPredictionLen() const {
        int len = 0;
        for (int i = 0; i < nrpparts.size(); ++i) {
            len += nrpparts[i].getAminoacidsPrediction().size();
        }

        return len;
    }
}
