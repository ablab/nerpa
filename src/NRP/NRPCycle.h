#ifndef NRPSMATCHER_NRPCYCLE_H
#define NRPSMATCHER_NRPCYCLE_H

#include "NRP.h"
#include "NRPLine.h"

namespace nrp {
    class NRPCycle : public NRP {

    public:
        NRPCycle(const std::string &file_name, const std::vector<std::string> &strformula,
                 const std::vector<aminoacid::Aminoacids::Aminoacid> &aminoacids, const std::vector<int> &position,
                 const std::string &graph, const std::string& extra_info);

        std::vector<Segment> containNRPsPart(nrpsprediction::NRPsPart predict_part) const override;

        NRPType getType() const override;

        std::vector<NRPLine> getLines() const override;
    };
}


#endif //NRPSMATCHER_NRPCYCLE_H
