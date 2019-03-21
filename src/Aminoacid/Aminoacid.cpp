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
    }

    bool Aminoacid::operator==(const Aminoacid &b) const {
        return formula == b.formula;
    }

    std::string Aminoacid::get_name() const {
        return  AminoacidInfo::AMINOACID_NAMES[aa];
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
}