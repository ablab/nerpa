#ifndef NRPSMATCHER_AMINOACIDPREDICTION_H
#define NRPSMATCHER_AMINOACIDPREDICTION_H

#include <vector>
#include <string>
#include "Aminoacid/Aminoacid.h"

namespace nrpsprediction {
    class AAdomainPrediction {
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
        bool is_repeatable=false;
        int pos;
        std::vector<AminoacidProb> aminoacid_prediction;
    public:
        AAdomainPrediction(int pos, std::vector<AminoacidProb> aminoacid_prediction, bool is_repeatable=false): pos(pos),
                                                                                      aminoacid_prediction(std::move(
                                                                                               aminoacid_prediction)),
                                                                                               is_repeatable(is_repeatable) {}
        bool contain(aminoacid::Aminoacid aminoacid) const;
        AminoacidProb getAminoacid(aminoacid::Aminoacid aminoacid) const;
        std::pair<int, int> getAmnAcidPos(aminoacid::Aminoacid aminoacid) const;
        bool is_rep() const {
            return is_repeatable;
        }
        std::vector<AminoacidProb> getAAPrediction() const {
            return aminoacid_prediction;
        }
    };
}


#endif //NRPSMATCHER_AMINOACIDPREDICTION_H
