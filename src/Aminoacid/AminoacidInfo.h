//
// Created by olga on 21.03.19.
//

#ifndef NRPSMATCHER_AMINOACIDINFO_H
#define NRPSMATCHER_AMINOACIDINFO_H

#include "Formula.h"

namespace aminoacid {
    class AminoacidInfo {
    public:
        enum AminoacidId {
            trp, ser, gly, uda, thr, dhp, gln, dab, arg, lys, ala_d, phe, val, cha, dhpg, phg, his, aeo,
            bmt, hse, met, ala, tcl, sal, allothr, b_ala, dhb, ile, end, leu, gua, hty, glu, bht, hpg, apa, pro, tyr,
            hyv, asn, cit, vol, cys, asp, dht, ahp, orn, apc, abu, aad, pip, dpg, none, AMINOACID_CNT
        };
        //int AMINOACID_CNT;
        static const std::string AMINOACID_NAMES[AMINOACID_CNT];
        static const Formula FORMULS[AMINOACID_CNT];

        AminoacidInfo() = default;
    };
}


#endif //NRPSMATCHER_AMINOACIDINFO_H
