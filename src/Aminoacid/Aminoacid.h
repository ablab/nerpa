#ifndef NRPSMATCHER_AMINOACIDS_H
#define NRPSMATCHER_AMINOACIDS_H

#include <string>
#include "Formula.h"
#include "Modification.h"
#include <vector>

namespace aminoacid {
    class Aminoacid {
    private:
        friend class ModificationTest;
    public:
        enum AminoacidId {
            trp, ser, gly, uda, thr, dhp, gln, dab, arg, lys, ala_d, phe, val, cha, dhpg, phg, his, aeo,
            bmt, hse, met, ala, tcl, sal, allothr, b_ala, dhb, ile, end, leu, gua, hty, glu, bht, hpg, apa, pro, tyr,
            hyv, asn, cit, vol, cys, asp, dht, ahp, orn, apc, abu, aad, pip, dpg, none, AMINOACID_CNT
        };
        static const std::string AMINOACID_NAMES[AMINOACID_CNT];
        static const Formula FORMULS[AMINOACID_CNT];

        explicit Aminoacid(std::string aminoacid_name);
        explicit Aminoacid(Formula formula);
        explicit Aminoacid(AminoacidId aid);
        Aminoacid();

        bool operator == (const Aminoacid& b) const;
        Formula operator - (const Aminoacid& b) const;
        friend std::ostream& operator << (std::ostream& out, const Aminoacid& a) {
            out << a.get_name();
            for (int i = 0; i < a.modifications.size(); ++i) {
                out << "+" << Modification::NAMES[a.modifications[i].getId()];
            }
            return out;
        }
        std::string get_name() const;
        void addModification(Modification m);

    private:
        Formula formula;
    public:
        const Formula &getFormula() const;

    private:
        AminoacidId aa;
        std::vector<Modification> modifications;
    };
}


#endif //NRPSMATCHER_AMINOACIDS_H
