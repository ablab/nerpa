//
// Created by olga on 01.03.19.
//

#ifndef NRPSMATCHER_PRISMPREDICTIONBUILDER_H
#define NRPSMATCHER_PRISMPREDICTIONBUILDER_H

#include <NRPsPrediction/NRPsPrediction.h>
#include "PredictionBuilderBase.h"

namespace nrpsprediction {
    class PrismPredictionBuilder : public PredictionBuilderBase {
    private:
        static const double EPS;
        std::vector<NRPsPart> nrpparts;
    public:
        void read_file(std::string file_name) override;

        NRPsPrediction getPrediction() override;
    };
}


#endif //NRPSMATCHER_PRISMPREDICTIONBUILDER_H
