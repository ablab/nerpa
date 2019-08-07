#include <fstream>
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
        FORMULS.push_back(Formula());
        MODIFICATION_CNT += 1;
        for (int i = 0; i < COEFFICIENT.size(); ++i) {
            COEFFICIENT[i].push_back(1);
        }

        in.close();
    }

    void ModificationInfo::init_AAMod(std::string filename) {
        std::fstream in(filename);
        std::string buffer;
        getline(in, buffer);

        while (getline(in, buffer)) {
            std::stringstream ss(buffer);
            std::string AA, modification;
            double coef;
            ss >> AA >> modification >> coef;
            int id = AminoacidInfo::getIdByNameId(AA);
            int mod_id = ModificationInfo::getIdByNameId(modification);
            if (id < AminoacidInfo::AMINOACID_CNT && mod_id < ModificationInfo::MODIFICATION_CNT) {
                COEFFICIENT[id][mod_id] = coef;
            }

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