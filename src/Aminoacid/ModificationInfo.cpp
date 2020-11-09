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
            MODIFICATION_CNT += 1;
            COEFFICIENT.push_back({coef_TT, coef_TF, coef_FT, coef_FF});

        }

        NAMES.push_back("empty");
        COEFFICIENT.push_back({0, 0, 0, 0});
        FORMULS.push_back(Formula(""));

        NAMES.push_back("");
        COEFFICIENT.push_back({0, 0, 0, 0});
        FORMULS.push_back(Formula(""));
        MODIFICATION_CNT = NAMES.size();

        in.close();
    }

    int ModificationInfo::getIdByNameId(std::string name) {
        for (int i = 0; i < MODIFICATION_CNT; ++i) {
            if (NAMES[i] == name) {
                return i;
            }
        }
        return MODIFICATION_CNT;
    }
}