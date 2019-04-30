#ifndef NRPSMATCHER_NRPSPART_H
#define NRPSMATCHER_NRPSPART_H

#include "AminoacidPrediction.h"
#include <string>

namespace nrpsprediction {
    class NRPsPart {
    private:
        std::string file_name = "";
        std::string orf = "";
        std::vector<AminoacidPrediction> aminoacids;
    public:
        NRPsPart(std::string file_name, std::string orf_name, int num,  const AminoacidPrediction& prediction);
        NRPsPart(std::string file_name, std::string orf_name);
        NRPsPart();
        std::string get_orf_name() const;
        std::string get_file_name() const;
        void add_prediction(int num, const AminoacidPrediction& prediction);
        std::vector<AminoacidPrediction> getAminoacidsPrediction() const;
    };
}


#endif //NRPSMATCHER_NRPSPART_H
