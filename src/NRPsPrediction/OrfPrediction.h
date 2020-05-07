#ifndef NRPSMATCHER_NRPSPART_H
#define NRPSMATCHER_NRPSPART_H

#include "AAdomainPrediction.h"
#include <string>

namespace nrpsprediction {
    class OrfPrediction {
    private:
        bool is_repeatable = false;
        std::string file_name = "";
        std::string orf = "";
        std::vector<AAdomainPrediction> aminoacids;
    public:
        OrfPrediction(std::string file_name, std::string orf_name, int num, const AAdomainPrediction& prediction, bool is_repeatable=false);
        OrfPrediction(std::string file_name, std::string orf_name, bool is_repeatable=false);
        OrfPrediction();
        std::string get_orf_name() const;
        std::string get_file_name() const;
        void add_prediction(int num, const AAdomainPrediction& prediction);
        std::vector<AAdomainPrediction> getAAdomainPrediction() const;
        bool is_rep() const  {
            return is_repeatable;
        }
    };
}


#endif //NRPSMATCHER_NRPSPART_H
