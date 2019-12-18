#include <cassert>
#include <iostream>
#include "Aminoacid.h"
#include "Formula.h"

namespace aminoacid {
    Aminoacid::Aminoacid(std::string aminoacid_name) {
        int was = 0;
        for (int i = 0; i < AminoacidInfo::AMINOACID_CNT; ++i) {
            if ( AminoacidInfo::AMINOACID_NAMES[i] == aminoacid_name) {
                this->aa = i;
                was = 1;
            }
        }
        if (was == 0) {
            std::cerr << aminoacid_name << "\n";
            assert(false);
        }

        formula = AminoacidInfo::FORMULS[aa];
    }

    Aminoacid::Aminoacid(Formula formula) {
        this->formula = formula;
        this->aa =  AminoacidInfo::AMINOACID_CNT - 1;
        this->possible_name = findAAname();
    }

    bool Aminoacid::operator==(const Aminoacid &b) const {
        return formula == b.formula;
    }

    std::string Aminoacid::get_name() const {
        return  AminoacidInfo::AMINOACID_NAMES[aa];
    }

    std::string Aminoacid::get_possible_name() const {
        if (!possible_name.empty()) {
            return possible_name;
        }

        return findAAname();
    }

    Aminoacid::Aminoacid(int aid) {
        this->aa = aid;
        formula = AminoacidInfo::FORMULS[aa];
    }

    Aminoacid::Aminoacid() {
        this->aa =  AminoacidInfo::AMINOACID_CNT - 1;
        formula = Formula();
    }

    Formula Aminoacid::operator-(const Aminoacid &b) const {
        return formula - b.formula;
    }

    void Aminoacid::addModification(Modification m) {
        formula += m.getFormula();
        modifications.push_back(m);
    }

    const Formula &Aminoacid::getFormula() const {
        return formula;
    }

    int Aminoacid::get_id() const {
        return aa;
    }

    std::string Aminoacid::findAAname() const {
        for (int i = 0; i < AminoacidInfo::AMINOACID_CNT; ++i) {
            if (this->formula == AminoacidInfo::FORMULS[i]) {
                return AminoacidInfo::AMINOACID_NAMES[i];
            }
        }

        for (int i = 0; i < AminoacidInfo::AMINOACID_CNT; ++i) {
            for (int j = 0; j < ModificationInfo::MODIFICATION_CNT; ++j) {
                if (this->formula - AminoacidInfo::FORMULS[i] == ModificationInfo::FORMULS[j]) {
                    return AminoacidInfo::AMINOACID_NAMES[i] + "+" + ModificationInfo::NAMES[j];
                }
            }
        }

        return AminoacidInfo::AMINOACID_NAMES[aa];
    }
}