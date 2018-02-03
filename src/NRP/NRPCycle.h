#ifndef NRPSMATCHER_NRPCYCLE_H
#define NRPSMATCHER_NRPCYCLE_H

#include "NRP.h"

namespace nrp {
    class NRPCycle : public NRP {

    public:
        NRPCycle(const std::string &file_name, const std::vector<std::string> &strformula,
                 const std::vector<aminoacid::Aminoacids::Aminoacid> &aminoacids, const std::vector<int> &position,
                 const std::string &graph);

        Match isCover(nrpsprediction::NRPsPrediction nrPsPrediction) override;
    };
}


#endif //NRPSMATCHER_NRPCYCLE_H
