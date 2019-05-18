//
// Created by olga on 01.02.19.
//

#include"Modification.h"

namespace aminoacid {
    const std::string Modification::NAMES[Modification::MODIFICATION_CNT] = {"methylation", "dimethylation", "demethylation",
                                                                             "hydration", "hydroxylation", "formylation",
                                                                             "phosphotylation", "acetylation", ""};
    const Formula Modification::FORMULS[Modification::MODIFICATION_CNT] = {Formula("CH2"), Formula("C2H4"), Formula("C-1H-2"),
                                                                           Formula("H2O"), Formula("O"), Formula("CO"),
                                                                           Formula("PH2O3"), Formula()};

    Formula Modification::getFormula() const {
        return formula;
    }

    Modification::Modification(Modification::ModificationId id) {
        this->id = id;
        formula = FORMULS[id];
    }

    Modification::Modification() {

    }

    Modification::ModificationId Modification::getId() const {
        return id;
    }

    Modification::Modification(Formula formula) {
        for (int i = 0; i < MODIFICATION_CNT; ++i) {
            if (formula == FORMULS[i]) {
                id = static_cast<ModificationId>(i);
                this->formula = formula;
                return;
            }
        }

        id = MODIFICATION_CNT;
        this->formula = Formula();
    }
}
