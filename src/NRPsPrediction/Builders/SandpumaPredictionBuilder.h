//
// Created by olga on 05.03.19.
//

#ifndef NRPSMATCHER_SANDPUMAPREDICTIONBUILDER_H
#define NRPSMATCHER_SANDPUMAPREDICTIONBUILDER_H

#include "PredictionBuilderBase.h"

namespace nrpsprediction {
    class SandpumaPredictionBuilder : public PredictionBuilderBase {
    private:
        //std::vector<AminoacidPrediction::AminoacidProb> parse_predictions(std::ifstream& in);
    public:
        static const std::string AMINOACID_NAMES[aminoacid::Aminoacid::AMINOACID_CNT];

        //parse ctg1_minowa_nrpspredoutput.txt file
        //expected groupded by orfs and sorted by num in one group
        void read_file(std::string file_name) override;
    };
}


#endif //NRPSMATCHER_SANDPUMAPREDICTIONBUILDER_H
