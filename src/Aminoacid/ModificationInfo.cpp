#include <fstream>
#include <Logger/logger.hpp>
#include "ModificationInfo.h"
#include "AminoacidInfo.h"

namespace aminoacid {
    int ModificationInfo::MODIFICATION_CNT;
    std::vector<std::string> ModificationInfo::NAMES;
    std::vector<Formula> ModificationInfo::FORMULS;
    std::vector<std::vector<double>> ModificationInfo::COEFFICIENT;

    void ModificationInfo::init(std::string filename) {
        std::fstream in(filename);
        std::string buffer;
        getline(in, buffer);

        while (getline(in, buffer)) {
            std::stringstream ss(buffer);
            std::string modification, formula;
            // coef_{PRED}{NRP}
            double coef_TT, coef_TF, coef_FT, coef_FF;
            ss >> modification >> coef_TT >> coef_TF >> coef_FT >> coef_FF;
            NAMES.push_back(modification);
            FORMULS.emplace_back("");
            COEFFICIENT.push_back({coef_TT, coef_TF, coef_FT, coef_FF});
        }

        NAMES.emplace_back("unk");
        COEFFICIENT.push_back({0, 0, 0, 0});
        FORMULS.emplace_back("");
        MODIFICATION_CNT = NAMES.size();

        in.close();
    }

    int ModificationInfo::getIdByNameId(std::string name) {
        for (int i = 0; i < MODIFICATION_CNT; ++i) {
            if (NAMES[i] == name) {
                return i;
            }
        }
        return MODIFICATION_CNT - 1;
    }

    double ModificationInfo::getCoefficientById(size_t id, bool in_pred, bool in_nrp) {
        size_t coef = 3; //FF
        if (in_pred && in_nrp) { //TT
            coef = 0;
        } else if (in_pred) { //TF
            coef = 1;
        } else if (in_nrp) { //FT
            coef = 2;
        }
        return COEFFICIENT[id][coef];
    }
}