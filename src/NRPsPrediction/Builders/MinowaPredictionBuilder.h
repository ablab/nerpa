//
// Created by olga on 23.02.19.
//

#ifndef NRPSMATCHER_MINOWAPREDICTIONBUILDER_H
#define NRPSMATCHER_MINOWAPREDICTIONBUILDER_H

#include "PredictionBuilderBase.h"

namespace nrpsprediction {
    class MinowaPredictionBuilder : public PredictionBuilderBase {
    private:
        std::vector<AAdomain_Prediction::AminoacidProb> parse_predictions(std::ifstream& in);
    public:

        //parse ctg1_minowa_nrpspredoutput.txt file
        //expected groupded by orfs and sorted by num in one group
        void read_file(std::string file_name) override;
    };
}


#endif //NRPSMATCHER_MINOWAPREDICTIONBUILDER_H
