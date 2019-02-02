#ifndef NRPSMATCHER_AMINOACIDPREDICTION_H
#define NRPSMATCHER_AMINOACIDPREDICTION_H

#include <vector>
#include <string>
#include "Aminoacid/Aminoacid.h"

namespace nrpsprediction {
    class AminoacidPrediction {
    public:
        struct AminoacidProb {
            AminoacidProb() {}

            aminoacid::Aminoacid aminoacid;
            double prob;

            AminoacidProb(std::string aacid, double prob): aminoacid(aacid), prob(prob) {}

            AminoacidProb(aminoacid::Aminoacid aacid, double prob): aminoacid(aacid), prob(prob) {}
        };
        static const double EPS;
    private:
        int pos;
        std::vector<AminoacidProb> aminoacid_prediction;
    public:
        AminoacidPrediction(int pos, std::string predict_aminoacids);
        std::pair<std::string, double> parse_token(std::string token);
        bool contain(aminoacid::Aminoacid aminoacid) const;
        AminoacidProb getAminoacid(aminoacid::Aminoacid aminoacid) const;
        std::pair<int, int> getAmnAcidPos(aminoacid::Aminoacid aminoacid) const;
        std::vector<AminoacidProb> getAAPrediction() const {
            return aminoacid_prediction;
        }
    };
}


#endif //NRPSMATCHER_AMINOACIDPREDICTION_H
