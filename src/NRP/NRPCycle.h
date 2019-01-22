#ifndef NRPSMATCHER_NRPCYCLE_H
#define NRPSMATCHER_NRPCYCLE_H

#include "NRP.h"

namespace nrp {
    class NRPCycle : public NRP {

    public:
        NRPCycle(const std::string &file_name, const std::vector<std::string> &strformula,
                 const std::vector<aminoacid::Aminoacids::Aminoacid> &aminoacids, const std::vector<int> &position,
                 const std::string &graph, const std::string& extra_info);

        Match isCover(const nrpsprediction::NRPsPrediction& nrPsPrediction) const override;

        std::vector<Segment> containNRPsPart(nrpsprediction::NRPsPart predict_part) const override;

        NRPType getType() const override;

    private:
        Match updateMatch(const nrpsprediction::NRPsPrediction& nrPsPrediction, Match match, int bg) const;
    };
}


#endif //NRPSMATCHER_NRPCYCLE_H
