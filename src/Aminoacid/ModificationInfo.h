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

    private:
 //       static void get_modification_list(std::string filename);

    };
}


#endif //NRPSMATCHER_MODIFICATIONINFO_H
