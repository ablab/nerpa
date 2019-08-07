//
// Created by olga on 01.02.19.
//

#include"Modification.h"
#include "ModificationInfo.h"

namespace aminoacid {
    Formula Modification::getFormula() const {
        return formula;
    }

    Modification::Modification(int id) {
        this->id = id;
        formula = ModificationInfo::FORMULS[id];
    }

    Modification::Modification() {

    }

    int Modification::getId() const {
        return id;
    }

    Modification::Modification(Formula formula) {
        for (int i = 0; i < ModificationInfo::MODIFICATION_CNT; ++i) {
            if (formula == ModificationInfo::FORMULS[i]) {
                id = i;
                this->formula = formula;
                return;
            }
        }

        id = ModificationInfo::MODIFICATION_CNT;
        this->formula = Formula();
    }
}
