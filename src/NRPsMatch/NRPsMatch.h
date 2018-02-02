#ifndef NRPSMATCHER_NRPSMATCH_H
#define NRPSMATCHER_NRPSMATCH_H

#include <vector>
#include <NRP/NRP.h>

namespace nrpsmatch {
    class NRPsMatch {
    private:
        nrp::NRP nrp;
        std::vector<nrpsprediction::NRPsPart> nrpParts;
        std::vector<int> parts_id;
        std::vector<int> parts_pos;
    public:
        NRPsMatch() = default;

        NRPsMatch(nrp::NRP nrp, std::vector<nrpsprediction::NRPsPart> nrpParts): nrp(nrp), nrpParts(nrpParts) {
            parts_id.resize(nrp.getLen(), -1);
            parts_pos.resize(nrp.getLen(), -1);
        }

        void match(int pos, int part_id, int part_pos);

        void print(std::ofstream& out);

        int score();

        std::vector<std::pair<int, int> > getMatchs();
    };
}
#endif //NRPSMATCHER_NRPSMATCH_H
