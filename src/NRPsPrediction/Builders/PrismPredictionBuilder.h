//
// Created by olga on 01.03.19.
//

#ifndef NRPSMATCHER_PRISMPREDICTIONBUILDER_H
#define NRPSMATCHER_PRISMPREDICTIONBUILDER_H

#include <NRPsPrediction/NRPsPrediction.h>
#include "PredictionBuilderBase.h"
#include <Json/json.hpp>

namespace nrpsprediction {
    class PrismPredictionBuilder : public PredictionBuilderBase {
    private:
        std::vector<AminoacidPrediction::AminoacidProb> parse_predictions(nlohmann::json predictions);
    public:
        static const std::string AMINOACID_NAMES[aminoacid::Aminoacid::AMINOACID_CNT];

        void read_file(std::string file_name) override;
    };
}


#endif //NRPSMATCHER_PRISMPREDICTIONBUILDER_H
