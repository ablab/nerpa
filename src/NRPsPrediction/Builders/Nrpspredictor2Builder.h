//
// Created by olga on 29.01.19.
//

#ifndef NRPSMATCHER_NRPSPREDICTOR2BUILDER_H
#define NRPSMATCHER_NRPSPREDICTOR2BUILDER_H

#include "PredictionBuilderBase.h"

namespace nrpsprediction {
    class Nrpspredictor2Builder : public PredictionBuilderBase {
    private:
        std::vector<NRPsPart> nrpparts;

        //orf = <prefix>_<orfname>_A<num>
        //return (<orfname>, <num>)
        std::pair<std::string, int> get_orf_name_and_order(std::string orf);
    public:
        //parse ctg1_nrpspredictor2_codes.txt file
        //expected groupded by orfs and sorted by num in one group
        void read_file(std::string file_name) override;

        NRPsPrediction getPrediction() override;
    };
}

#endif //NRPSMATCHER_NRPSPREDICTOR2BUILDER_H
