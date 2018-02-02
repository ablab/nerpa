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
        NRPsPart(std::string file_name, std::string orf_name, int num,  std::string predict_aminoacids);
        NRPsPart(std::string file_name, std::string orf_name);
        std::string get_orf_name();
        std::string get_file_name();
        void add_prediction(int num, std::string predict_aminoacids);
        std::vector<AminoacidPrediction> getAminoacidsPrediction();
    };
}


#endif //NRPSMATCHER_NRPSPART_H
