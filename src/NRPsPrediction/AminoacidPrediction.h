#ifndef NRPSMATCHER_AMINOACIDPREDICTION_H
#define NRPSMATCHER_AMINOACIDPREDICTION_H

#include <vector>
#include <string>
#include "../Aminoacid/Aminoacids.h"

namespace nrpsprediction {
    class AminoacidPrediction {
    public:
        struct AminoacidProb {
            aminoacid::Aminoacids::Aminoacid aminoacid;
            double prob;

            AminoacidProb(std::string aacid, double prob) {
                this->prob = prob;
                aminoacid = aminoacid::Aminoacids::get_aminoacid(aacid);
            }

            AminoacidProb(aminoacid::Aminoacids::Aminoacid aacid, double prob) {
                this->prob = prob;
                aminoacid = aacid;
            }
        };
        static const double EPS;
    private:
        int pos;
        std::vector<AminoacidProb> aminoacid_prediction;
    public:
        AminoacidPrediction(int pos, std::vector<AminoacidProb> aminoacid_prediction): pos(pos),
                                                                                       aminoacid_prediction(aminoacid_prediction) {}
        AminoacidPrediction(int pos, std::string predict_aminoacids);
        std::pair<std::string, double> parse_token(std::string token);
        bool contain(aminoacid::Aminoacids::Aminoacid aminoacid);
        double getScore(aminoacid::Aminoacids::Aminoacid aminoacid);
        AminoacidProb getAminoacid(aminoacid::Aminoacids::Aminoacid aminoacid);
        std::pair<int, int> getAmnAcidPos(aminoacid::Aminoacids::Aminoacid aminoacid);
        std::vector<AminoacidProb> getAAPrediction() {
            return aminoacid_prediction;
        }
    };
}


#endif //NRPSMATCHER_AMINOACIDPREDICTION_H
