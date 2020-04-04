#ifndef NRPSMATCHER_NRPSPREDICTION_H
#define NRPSMATCHER_NRPSPREDICTION_H

#include "ORF_Prediction.h"
#include <iostream>

namespace nrpsprediction {
    //structure for store NRPs predictions
    class BGC_Prediction {
    private:
        std::vector<ORF_Prediction> nrpparts;
        std::vector<ORF_Prediction> short_parts;
    public:
        BGC_Prediction() = default;
        BGC_Prediction(std::vector<ORF_Prediction> nrpparts): nrpparts(nrpparts) {}
        BGC_Prediction(std::vector<ORF_Prediction> nrpparts, std::vector<ORF_Prediction> short_parts):
                nrpparts(nrpparts), short_parts(short_parts) {}
        const std::vector<ORF_Prediction>& getNrpsParts() const;
        const std::vector<ORF_Prediction>& getShortParts() const;

        int getSumPredictionLen() const;
    };
}


#endif //NRPSMATCHER_NRPSPREDICTION_H
