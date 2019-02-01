//
// Created by olga on 01.02.19.
//

#ifndef NRPSMATCHER_MODIFICATION_H
#define NRPSMATCHER_MODIFICATION_H

#include "Formula.h"

namespace aminoacid {
    class Modification {
    public:
        enum ModificationId {OH, me3, MODIFICATION_CNT};
        static const std::string NAMES[MODIFICATION_CNT];
        static const Formula FORMULS[MODIFICATION_CNT];

    private:
        Formula formula;
        ModificationId id;
    };
}

#endif //NRPSMATCHER_MODIFICATION_H
