#ifndef NRPSMATCHER_NRPLINE_H
#define NRPSMATCHER_NRPLINE_H

#include "NRP.h"

namespace nrp {
    class NRPLine : public NRP {
    public:
        NRPLine(const std::string &file_name, const std::vector<std::string> &strformula,
                const std::vector<aminoacid::Aminoacids::Aminoacid> &aminoacids, const std::vector<int> &position,
                const std::string &graph);

        Match isCover(nrpsprediction::NRPsPrediction nrPsPrediction) override;

        std::vector<Segment> containNRPsPart(nrpsprediction::NRPsPart predict_part) override;
    };
}


#endif //NRPSMATCHER_NRPLINE_H
