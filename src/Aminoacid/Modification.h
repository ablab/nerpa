//
// Created by olga on 01.02.19.
//

#ifndef NRPSMATCHER_MODIFICATION_H
#define NRPSMATCHER_MODIFICATION_H

#include "Formula.h"

namespace aminoacid {
    class Modification {
    public:
        enum ModificationId {methylation, dimethylation, demethylation,
            hydration, hydroxylation, formylation,
            phosphotylation, acetylation,
            dedimethylation,
            dehydration, dehydroxylation, deformylation,
            dephosphotylation, deacetylation,
            empty, MODIFICATION_CNT};
        static const std::string NAMES[MODIFICATION_CNT];
        static const Formula FORMULS[MODIFICATION_CNT];

        Formula getFormula() const;
        ModificationId getId() const;

        explicit Modification(ModificationId id);
        explicit Modification(Formula formula);
        Modification();
    private:
        Formula formula;
        ModificationId id;
    };
}

#endif //NRPSMATCHER_MODIFICATION_H
