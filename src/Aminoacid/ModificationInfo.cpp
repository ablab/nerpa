#include <fstream>
#include "ModificationInfo.h"

namespace aminoacid {
    int ModificationInfo::MODIFICATION_CNT;
    std::vector<std::string> ModificationInfo::NAMES;
    std::vector<Formula> ModificationInfo::FORMULS;
    std::vector<std::vector<double>> ModificationInfo::COEFFICIENT;

    void ModificationInfo::init(std::string filename) {

    }
}