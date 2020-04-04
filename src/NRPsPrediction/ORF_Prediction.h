#ifndef NRPSMATCHER_NRPSPART_H
#define NRPSMATCHER_NRPSPART_H

#include "AAdomain_Prediction.h"
#include <string>

namespace nrpsprediction {
    class ORF_Prediction {
    private:
        std::string file_name = "";
        std::string orf = "";
        std::vector<AAdomain_Prediction> aminoacids;
    public:
        ORF_Prediction(std::string file_name, std::string orf_name, int num, const AAdomain_Prediction& prediction);
        ORF_Prediction(std::string file_name, std::string orf_name);
        ORF_Prediction();
        std::string get_orf_name() const;
        std::string get_file_name() const;
        void add_prediction(int num, const AAdomain_Prediction& prediction);
        std::vector<AAdomain_Prediction> getAminoacidsPrediction() const;
    };
}


#endif //NRPSMATCHER_NRPSPART_H
