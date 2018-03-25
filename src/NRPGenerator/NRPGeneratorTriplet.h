#ifndef NRPSMATCHER_NRPGENERATORTRIPLET_H
#define NRPSMATCHER_NRPGENERATORTRIPLET_H

#include "NRPGenerator.h"

namespace nrp_generator {
typedef aminoacid::Aminoacids::Aminoacid Aminoacid;
    class NRPGeneratorTriplet : public NRPGenerator {
    private:
        std::vector<std::vector<int> > graph;
        std::vector<int> cntVert;
        std::pair<Aminoacid, Aminoacid> getVertexAA(int id);
        int getVertexId(Aminoacid aa1, Aminoacid aa2);

        nrp::NRP* generateNRPCycle(int len);
        nrp::NRP* generateNRPBranchCycle(int len);
        nrp::NRP* generateNRPLine(int len);


        Aminoacid genAA(int& v);
    public:
        explicit NRPGeneratorTriplet(const std::vector<nrp::NRP *> &mols);
        nrp::NRP *generate(nrp::NRP::NRPType type, int len) override;

        int genStartVert();
    };
}


#endif //NRPSMATCHER_NRPGENERATORTRIPLET_H
