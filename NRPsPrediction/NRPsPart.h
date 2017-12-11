#ifndef NRPSMATCHER_NRPSPART_H
#define NRPSMATCHER_NRPSPART_H

#include "AminoacidPrediction.h"
#include <string>

namespace nrpsprediction {
    class NRPsPart {
    private:
        std::string file_name;
        std::string orf;
        std::vector<AminoacidPrediction> aminoacids;
    public:
        std::string get_orf_name();
        void add_prediction(int num, std::string predict_aminoacids)
    };
}


#endif //NRPSMATCHER_NRPSPART_H
