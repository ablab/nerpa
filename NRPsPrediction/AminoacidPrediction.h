#ifndef NRPSMATCHER_AMINOACIDPREDICTION_H
#define NRPSMATCHER_AMINOACIDPREDICTION_H

#include <vector>
#include <string>
#include "../Aminoacid/Aminoacids.h"

namespace nrpsprediction {
    class AminoacidPrediction {
    private:
        struct AminoacidProb {
            aminoacid::Aminoacids::Aminoacid aminoacid;
            double prob;

            AminoacidProb(std::string aacid, double prob) {
                this->prob = prob;
                aminoacid = aminoacid::Aminoacids::getAminoacid(aacid);
            }
        };
        const double EPS = 1e-4;
    private:
        int pos;
        std::vector<AminoacidProb> aminoacid_prediction;
    public:
        AminoacidPrediction(int pos, std::string predict_aminoacids);

        std::pair<std::string, double> parse_token(std::string token);
    };
}


#endif //NRPSMATCHER_AMINOACIDPREDICTION_H
