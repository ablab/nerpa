#ifndef NRPSMATCHER_NRPGENERATOR_H
#define NRPSMATCHER_NRPGENERATOR_H

#include "NRP/NRP.h"

namespace nrp_generator {
    class NRPGenerator {
    private:
        int sumAA = 0;
        std::vector<int> cntAA;

        nrp::NRP* generateNRPCycle(int len);
        nrp::NRP* generateNRPBranchCycle(int len);
        nrp::NRP* generateNRPLine(int len);
    public:
        explicit NRPGenerator(std::vector<nrp::NRP *> mols);

        virtual nrp::NRP* generate(nrp::NRP::NRPType type, int len);

        aminoacid::Aminoacid::AminoacidId genAA();
    };
}


#endif //NRPSMATCHER_NRPGENERATOR_H
