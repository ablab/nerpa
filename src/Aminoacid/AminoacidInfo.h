//
// Created by olga on 21.03.19.
//

#ifndef NRPSMATCHER_AMINOACIDINFO_H
#define NRPSMATCHER_AMINOACIDINFO_H

#include "Formula.h"
#include <vector>

namespace aminoacid {
    class AminoacidInfo {
    public:
        static double DEFAULT_LOG_P;
        static int AMINOACID_CNT;
        static std::vector<std::string> AMINOACID_NAMES;
        static std::vector<std::string> NAME_ID;
        static std::vector<Formula> FORMULS;
        static std::vector<double> LogP;

        static void init(std::string filename, std::string predictor, double default_logp);
        static int getIdByNameId(std::string name);
    };
}


#endif //NRPSMATCHER_AMINOACIDINFO_H
