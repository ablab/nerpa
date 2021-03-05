#ifndef NRPSMATCHER_MODIFICATIONINFO_H
#define NRPSMATCHER_MODIFICATIONINFO_H

#include "Formula.h"

namespace aminoacid {
    class ModificationInfo {
    public:
        static int MODIFICATION_CNT;
        static std::vector<std::string> NAMES;
        static std::vector<Formula> FORMULS; // TODO: Do we still need this?
        static std::vector<std::vector<double>> COEFFICIENT;

        static void init(std::string filename);
//        static void init_AAMod(std::string filename);
        static int getIdByNameId(std::string name);
        static double getCoefficientById(size_t id, bool in_pred, bool in_nrp);
    };
}


#endif //NRPSMATCHER_MODIFICATIONINFO_H
