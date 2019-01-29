//
// Created by olga on 29.01.19.
//

#ifndef NRPSMATCHER_PREDICTIONBUILDERBASE_H
#define NRPSMATCHER_PREDICTIONBUILDERBASE_H

#include <NRPsPrediction/NRPsPrediction.h>

namespace nrpsprediction {
    class PredictionBuilderBase {
    public:
        virtual void read_file(std::string file_name) = 0;

        virtual NRPsPrediction getPrediction() = 0;
    };
}


#endif //NRPSMATCHER_PREDICTIONBUILDERBASE_H
