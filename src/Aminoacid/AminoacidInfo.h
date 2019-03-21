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
        static int AMINOACID_CNT;
        static std::vector<std::string> AMINOACID_NAMES;
        static std::vector<Formula> FORMULS;

        static void init(std::string filename, std::string predictor);
    };
}


#endif //NRPSMATCHER_AMINOACIDINFO_H
