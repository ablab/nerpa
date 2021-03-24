//
// Created by tag on 21/03/2020.
//

#ifndef NERPA_MONOMERINFO_H
#define NERPA_MONOMERINFO_H

#include <vector>
#include <string>
#include <unordered_map>
#include "Aminoacid.h"

namespace aminoacid {
    class MonomerInfo {
    public:
        static double DEFAULT_LOG_P;
        static std::unordered_map<std::string, int> MONOMER_TO_AA;
        static std::unordered_map<std::string, std::vector<int>> MONOMER_TO_MODIFICATIONS;
        static std::unordered_map<int, double> LogP;

        static void init(const std::string& filename, const std::string& logp_filename, double default_logp);
        static double getLogP(std::string& code);
        static double getLogP(int AAid);
        static int getAAIdByCode(std::string& code);
        static Aminoacid getAAByCode(std::string& code);
    };
}

#endif //NERPA_MONOMERINFO_H
