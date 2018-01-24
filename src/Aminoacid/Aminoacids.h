#ifndef NRPSMATCHER_AMINOACIDS_H
#define NRPSMATCHER_AMINOACIDS_H

#include <string>

namespace aminoacid {
    class Aminoacids {
    public:
        enum Aminoacid {trp, ser, gly, uda, thr, dhp, gln, dab, arg, lys, ala_d, phe, val, cha, dhpg, phg, his, aeo,
            bmt, hse, met, ala, tcl, sal, allothr, b_ala, dhb, ile, end, leu, gua, hty, glu, bht, hpg, apa, pro, tyr,
            hyv, asn, cit, vol, cys, asp, dht, ahp, orn, apc, abu, aad, pip, dpg, AMINOACID_CNT};

        static const std::string AMINOACID_NAMES[AMINOACID_CNT];

        static  const std::string FORMULS[AMINOACID_CNT];

        static Aminoacid get_aminoacid(std::string aminoacid_name);
        static Aminoacid get_aminoacid_from_formula(std::string fotmula);
        static bool same(Aminoacid a, Aminoacid b);
    };

    }


#endif //NRPSMATCHER_AMINOACIDS_H
