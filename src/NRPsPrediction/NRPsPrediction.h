#ifndef NRPSMATCHER_NRPSPREDICTION_H
#define NRPSMATCHER_NRPSPREDICTION_H

#include "NRP/NRP.h"
#include "NRPsMatch/NRPsMatch.h"
#include "NRPsPart.h"

namespace nrpsprediction {
    using namespace nrpsmatch;
    //structure for store NRPs predictions
    class NRPsPrediction {
    private:
        std::vector<NRPsPart> nrpparts;

        //orf = <prefix>_<orfname>_A<num>
        //return (<orfname>, <num>)
        std::pair<std::string, int> get_orf_name_and_order(std::string orf);


        struct Segment {
            int l;
            int r;
            int part_id;
            bool rev;

            Segment(){}
            Segment(int l, int r, int id, bool rev): l(l), r(r), part_id(id), rev(rev) {}

            bool operator < (Segment b) {
                return l < b.l;
            }
        };

        nrpsmatch::NRPsMatch isCoverLine(std::vector<Segment>& segments, nrp::NRP nrp,
                         const std::vector<int>& toSmallId, const std::vector<int>& toBigId);

    public:
        //parse ctg1_nrpspredictor2_codes.txt file
        //expected groupded by orfs and sorted by num in one group
        void read_file(std::string file_name);
        std::vector<NRPsPart> getNrpsParts();

        nrpsmatch::NRPsMatch isCover(nrp::NRP nrp);
    };
}


#endif //NRPSMATCHER_NRPSPREDICTION_H
