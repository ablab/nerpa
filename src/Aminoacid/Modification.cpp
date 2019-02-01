//
// Created by olga on 01.02.19.
//

#include"Modification.h"

namespace aminoacid {
    const std::string Modification::NAMES[Modification::MODIFICATION_CNT] = {"OH", "3me"};
    const Formula Modification::FORMULS[Modification::MODIFICATION_CNT] = {Formula("O-2"), Formula("CH2")};
}
