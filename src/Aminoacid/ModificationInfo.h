#ifndef NRPSMATCHER_MODIFICATIONINFO_H
#define NRPSMATCHER_MODIFICATIONINFO_H

#include "Formula.h"

namespace aminoacid {
    class ModificationInfo {
    public:
        static int MODIFICATION_CNT;
        static std::vector<std::string> NAMES;
        static std::vector<Formula> FORMULS;
        static std::vector<std::vector<double>> COEFFICIENT;

        static void init(std::string filename);
        static void init_AAMod(std::string filename);
        static int getIdByNameId(std::string name);
    };
}


#endif //NRPSMATCHER_MODIFICATIONINFO_H
