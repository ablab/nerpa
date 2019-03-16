//
// Created by olga on 29.01.19.
//

#ifndef NRPSMATCHER_NRPSPREDICTOR2BUILDER_H
#define NRPSMATCHER_NRPSPREDICTOR2BUILDER_H

#include "PredictionBuilderBase.h"

namespace nrpsprediction {
    class Nrpspredictor2Builder : public PredictionBuilderBase {
    private:
        static std::pair<std::string, double> parse_token(std::string token);
    public:
        static std::vector<AminoacidPrediction::AminoacidProb> parse_predictions(std::string predictions);

        //parse ctg1_nrpspredictor2_codes.txt file
        //expected groupded by orfs and sorted by num in one group
        void read_file(std::string file_name) override;
    };
}

#endif //NRPSMATCHER_NRPSPREDICTOR2BUILDER_H
