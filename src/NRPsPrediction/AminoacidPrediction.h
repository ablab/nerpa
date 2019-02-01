#ifndef NRPSMATCHER_AMINOACIDPREDICTION_H
#define NRPSMATCHER_AMINOACIDPREDICTION_H

#include <vector>
#include <string>
#include "Aminoacid/Aminoacid.h"

namespace nrpsprediction {
    class AminoacidPrediction {
    public:
        struct AminoacidProb {
            aminoacid::Aminoacid::AminoacidId aminoacid;
            double prob;

            AminoacidProb(std::string aacid, double prob) {
                this->prob = prob;
                aminoacid = aminoacid::Aminoacid::get_aminoacid(aacid);
            }

            AminoacidProb(aminoacid::Aminoacid::AminoacidId aacid, double prob) {
                this->prob = prob;
                aminoacid = aacid;
            }
        };
        static const double EPS;
    private:
        int pos;
        std::vector<AminoacidProb> aminoacid_prediction;
    public:
        AminoacidPrediction(int pos, std::string predict_aminoacids);
        std::pair<std::string, double> parse_token(std::string token);
        bool contain(aminoacid::Aminoacid::AminoacidId aminoacid);
        AminoacidProb getAminoacid(aminoacid::Aminoacid::AminoacidId aminoacid) const;
        std::pair<int, int> getAmnAcidPos(aminoacid::Aminoacid::AminoacidId aminoacid) const;
        std::vector<AminoacidProb> getAAPrediction() {
            return aminoacid_prediction;
        }
    };
}


#endif //NRPSMATCHER_AMINOACIDPREDICTION_H
