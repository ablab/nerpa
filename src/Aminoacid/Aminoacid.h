#ifndef NRPSMATCHER_AMINOACIDS_H
#define NRPSMATCHER_AMINOACIDS_H

#include <string>
#include "Formula.h"

namespace aminoacid {
    class Aminoacid {
    public:
        enum AminoacidId {
            trp, ser, gly, uda, thr, dhp, gln, dab, arg, lys, ala_d, phe, val, cha, dhpg, phg, his, aeo,
            bmt, hse, met, ala, tcl, sal, allothr, b_ala, dhb, ile, end, leu, gua, hty, glu, bht, hpg, apa, pro, tyr,
            hyv, asn, cit, vol, cys, asp, dht, ahp, orn, apc, abu, aad, pip, dpg, me3_glu, OH_asn, none, AMINOACID_CNT
        };
        static const std::string AMINOACID_NAMES[AMINOACID_CNT];
        static const Formula FORMULS[AMINOACID_CNT];

        explicit Aminoacid(std::string aminoacid_name);
        explicit Aminoacid(Formula formula);
        explicit Aminoacid(AminoacidId aid);
        Aminoacid();

        static AminoacidId get_aminoacid(std::string aminoacid_name);

        static AminoacidId get_aminoacid_from_formula(std::string fotmula);

        static bool same(AminoacidId a, AminoacidId b);
        bool operator == (const Aminoacid& b) const;
        std::string get_name() const;

    private:
        Formula formula;
        AminoacidId aa;
    };
}


#endif //NRPSMATCHER_AMINOACIDS_H
