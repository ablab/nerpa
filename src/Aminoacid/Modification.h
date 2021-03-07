//
// Created by olga on 01.02.19.
//

#ifndef NRPSMATCHER_MODIFICATION_H
#define NRPSMATCHER_MODIFICATION_H

#include "Formula.h"

namespace aminoacid {
    class Modification {
    public:
        Formula getFormula() const;
        int getId() const;

        explicit Modification(int id);
        explicit Modification(Formula formula);
        Modification();

//        double getScore(int AAid);
    private:
        Formula formula;
        int id;
    };
}

#endif //NRPSMATCHER_MODIFICATION_H
