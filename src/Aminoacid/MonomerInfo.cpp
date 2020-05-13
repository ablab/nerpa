//
// Created by tag on 21/03/2020.
//

#include <fstream>
#include "MonomerInfo.h"
#include "AminoacidInfo.h"
#include "ModificationInfo.h"

namespace aminoacid {

    std::unordered_map<std::string, int> MonomerInfo::MONOMER_TO_AA;
    std::unordered_map<std::string, std::vector<int>> MonomerInfo::MONOMER_TO_MODIFICATIONS;

    void MonomerInfo::init(std::string& filename) {
        std::fstream in(filename);
        std::string buffer;
        getline(in, buffer); // skip header

        while (getline(in, buffer)) {
            std::stringstream ss(buffer);
            std::string code, nameId, modifications;
            ss >> code >> nameId >> modifications;
            MONOMER_TO_AA[code] = aminoacid::AminoacidInfo::getIdByNameId(nameId);
            std::vector<int> modificationIds;
            if (modifications != "-") {
                std::stringstream ss_(modifications);
                std::string single_mod;
                while (std::getline(ss_, single_mod, '+')) {
                    modificationIds.push_back(aminoacid::ModificationInfo::getIdByNameId(single_mod));
                }
            }
            MONOMER_TO_MODIFICATIONS[code] = modificationIds;
        }
    }

    int MonomerInfo::getAAIdByCode(std::string& code) {
        if (MONOMER_TO_AA.find(code) != MONOMER_TO_AA.end()) {
            return MONOMER_TO_AA[code];
        }
        return AminoacidInfo::AMINOACID_CNT - 1;
    }

    Aminoacid MonomerInfo::getAAByCode(std::string &code) {
        int aid = getAAIdByCode(code);
        Aminoacid res = Aminoacid(aid);
        if (MONOMER_TO_MODIFICATIONS.find(code) != MONOMER_TO_MODIFICATIONS.end()) {
            for (int mid : MONOMER_TO_MODIFICATIONS[code]) {
                Modification mod(mid);
                res.addModification(mod);
            }
        }
        return res;
    }
}