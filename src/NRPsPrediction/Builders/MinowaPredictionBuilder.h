//
// Created by olga on 23.02.19.
//

#ifndef NRPSMATCHER_MINOWAPREDICTIONBUILDER_H
#define NRPSMATCHER_MINOWAPREDICTIONBUILDER_H

#include "PredictionBuilderBase.h"

namespace nrpsprediction {
    class MinowaPredictionBuilder : public PredictionBuilderBase {
    private:
        static const double EPS;
        std::vector<NRPsPart> nrpparts;
        std::vector<AminoacidPrediction::AminoacidProb> parse_predictions(std::ifstream& in);
        aminoacid::Aminoacid::AminoacidId getAAbyName(std::string s);
    public:
        static const std::string AMINOACID_NAMES[aminoacid::Aminoacid::AMINOACID_CNT];

        //parse ctg1_minowa_nrpspredoutput.txt file
        //expected groupded by orfs and sorted by num in one group
        void read_file(std::string file_name) override;

        NRPsPrediction getPrediction() override;
    };
}


#endif //NRPSMATCHER_MINOWAPREDICTIONBUILDER_H
