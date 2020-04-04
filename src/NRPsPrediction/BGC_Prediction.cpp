#include <fstream>
#include <sstream>
#include <iostream>
#include <algorithm>
#include <Logger/logger.hpp>
#include "BGC_Prediction.h"

namespace nrpsprediction {
    const std::vector<ORF_Prediction>& BGC_Prediction::getNrpsParts() const {
        return nrpparts;
    }

    const std::vector<ORF_Prediction>& BGC_Prediction::getShortParts() const {
        return short_parts;
    }

    int BGC_Prediction::getSumPredictionLen() const {
        int len = 0;
        for (int i = 0; i < nrpparts.size(); ++i) {
            len += nrpparts[i].getAminoacidsPrediction().size();
        }

        return len;
    }
}
