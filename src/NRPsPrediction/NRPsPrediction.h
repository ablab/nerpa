#ifndef NRPSMATCHER_NRPSPREDICTION_H
#define NRPSMATCHER_NRPSPREDICTION_H

#include "NRPsPart.h"
#include <iostream>

namespace nrpsprediction {
    //structure for store NRPs predictions
    class NRPsPrediction {
    private:
        std::vector<NRPsPart> nrpparts;
        std::vector<NRPsPart> short_parts;
    public:
        NRPsPrediction() = default;
        NRPsPrediction(std::vector<NRPsPart> nrpparts): nrpparts(nrpparts) {}
        NRPsPrediction(std::vector<NRPsPart> nrpparts, std::vector<NRPsPart> short_parts):
                nrpparts(nrpparts), short_parts(short_parts) {}
        const std::vector<NRPsPart>& getNrpsParts() const;
        const std::vector<NRPsPart>& getShortParts() const;

        int getSumPredictionLen() const;
    };
}


#endif //NRPSMATCHER_NRPSPREDICTION_H
