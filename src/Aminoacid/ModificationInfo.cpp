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
        COEFFICIENT.resize(AminoacidInfo::AMINOACID_CNT);
        std::fstream in(filename);
        std::string buffer;
        getline(in, buffer);

        while (getline(in, buffer)) {
            std::stringstream ss(buffer);
            std::string modification, formula;
            double coef;
            ss >> modification >> formula >> coef;
            NAMES.push_back(modification);
            FORMULS.push_back(Formula(formula));
            MODIFICATION_CNT += 1;
            for (int i = 0; i < COEFFICIENT.size(); ++i) {
                COEFFICIENT[i].push_back(coef);
            }
        }

        NAMES.push_back("empty");
        NAMES.push_back("unk");
        FORMULS.push_back(Formula());
        FORMULS.push_back(Formula());
        MODIFICATION_CNT += 1;
        for (int i = 0; i < COEFFICIENT.size(); ++i) {
            COEFFICIENT[i].push_back(1);
            COEFFICIENT[i].push_back(-1);
        }

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