#ifndef NRPSMATCHER_NRPSPREDICTION_H
#define NRPSMATCHER_NRPSPREDICTION_H

#include "NRPsPart.h"
#include <iostream>

namespace nrpsprediction {
    //structure for store NRPs predictions
    class NRPsPrediction {
    private:
        std::vector<NRPsPart> nrpparts;
    public:
        NRPsPrediction() = default;
        NRPsPrediction(std::vector<NRPsPart> nrpparts): nrpparts(nrpparts) {}
        std::vector<NRPsPart> getNrpsParts() const;
    };
}


#endif //NRPSMATCHER_NRPSPREDICTION_H
