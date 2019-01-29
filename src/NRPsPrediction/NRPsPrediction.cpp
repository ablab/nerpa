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
}
