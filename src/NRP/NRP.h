#ifndef NRPSMATCHER_NRP_H
#define NRPSMATCHER_NRP_H

#include <vector>
#include "../Aminoacid/Aminoacids.h"
#include "../NRPsPrediction/NRPsPart.h"

namespace nrp {
    class NRP {
    private:
        friend class ContainNRPsTest;
        std::vector <aminoacid::Aminoacids::Aminoacid> aminoacids;
    public:
        void parse_fragment_graph(std::string fragment_graph);
        bool containNRPsPart(nrpsprediction::NRPsPart predict_part);
        void print();
    };
}


#endif //NRPSMATCHER_NRP_H
