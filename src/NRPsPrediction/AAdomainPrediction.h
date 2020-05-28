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
            std::vector<aminoacid::Modification> modificatins;
            double prob;
            AminoacidProb(std::string aacid, double prob): aminoacid(aacid), prob(prob) {}
            AminoacidProb(aminoacid::Aminoacid aacid, double prob): aminoacid(aacid), prob(prob) {}
        };
        static const double EPS;
    private:
        bool is_repeatable=false;
        int pos;
        std::vector<AminoacidProb> aminoacid_prediction;
        std::vector<aminoacid::Modification> modificatins;
    public:
        AAdomainPrediction(int pos, std::vector<AminoacidProb> aminoacid_prediction, bool is_repeatable=false, std::vector<aminoacid::Modification> mods={}, aminoacid::Aminoacid::Configuation configuation=aminoacid::Aminoacid::NA): pos(pos),
                                                                                      aminoacid_prediction(std::move(
                                                                                               aminoacid_prediction)),
                                                                                               is_repeatable(is_repeatable),
                                                                                               modificatins(std::move(mods)) {
            for (int i = 0; i < this->modificatins.size(); ++i) {
                for (int j = 0; j < this->aminoacid_prediction.size(); ++j) {
                    this->aminoacid_prediction[j].modificatins.push_back(this->modificatins[i]);
                }
            }

            for (int j = 0; j < this->aminoacid_prediction.size(); ++j) {
                this->aminoacid_prediction[j].aminoacid.setConfiguration(configuation);
            }
        }

        bool contain(aminoacid::Aminoacid aminoacid) const;
        AminoacidProb getAminoacid(aminoacid::Aminoacid aminoacid) const;
        std::pair<int, int> getAmnAcidPos(aminoacid::Aminoacid aminoacid) const;
        bool is_rep() const {
            return is_repeatable;
        }
        std::vector<AminoacidProb> getAAPrediction() const {
            return aminoacid_prediction;
        }

        std::vector<aminoacid::Modification> get_modifications() const {
            return modificatins;
        }
    };
}


#endif //NRPSMATCHER_AMINOACIDPREDICTION_H
