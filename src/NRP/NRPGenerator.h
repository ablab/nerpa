#ifndef NRPSMATCHER_NRPGENERATOR_H
#define NRPSMATCHER_NRPGENERATOR_H


#include "NRP.h"

namespace nrp {
    class NRPGenerator {
    private:
        int sumAA = 0;
        std::vector<int> cntAA;

        nrp::NRP* generateNRPCycle(int len);
        nrp::NRP* generateNRPBranchCycle(int len);
        nrp::NRP* generateNRPLine(int len);
    public:
        NRPGenerator(std::vector<nrp::NRP *> mols);

        nrp::NRP* generate(nrp::NRP::NRPType type, int len);

        aminoacid::Aminoacids::Aminoacid genAA();
    };
}


#endif //NRPSMATCHER_NRPGENERATOR_H
