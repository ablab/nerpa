//
// Created by olga on 29.01.19.
//

#ifndef NRPSMATCHER_PREDICTIONBUILDERBASE_H
#define NRPSMATCHER_PREDICTIONBUILDERBASE_H

#include <NRPsPrediction/BgcPrediction.h>

namespace nrpsprediction {
    class PredictionBuilderBase {
    protected:
        static const double EPS;
        std::vector<OrfPrediction> nrpparts;
        std::vector<OrfPrediction> short_parts;
        static int getAAbyName(std::string s);

        //orf = <prefix>_<orfname>_A<num>
        //return (<orfname>, <num>)
        std::pair<std::string, int> get_orf_name_and_order(std::string orf);
    public:
        virtual void read_file(std::string file_name) = 0;

        virtual BgcPrediction getPrediction();
    };
}


#endif //NRPSMATCHER_PREDICTIONBUILDERBASE_H
