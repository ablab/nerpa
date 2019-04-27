#include <fstream>
#include <sstream>
#include <iostream>
#include <algorithm>
#include <Logger/logger.hpp>
#include "NRPsPrediction.h"

namespace nrpsprediction {
    std::vector<NRPsPart> NRPsPrediction::getNrpsParts() const {
        return nrpparts;
    }

    int NRPsPrediction::getSumPredictionLen() const {
        int len = 0;
        for (int i = 0; i < nrpparts.size(); ++i) {
            len += nrpparts[i].getAminoacidsPrediction().size();
        }

        return len;
    }
}
