#ifndef NRPSMATCHER_NRPLINE_H
#define NRPSMATCHER_NRPLINE_H

#include "NRP.h"

namespace nrp {
    class NRPLine : public NRP {
    public:
        NRPLine(const std::string &file_name, const std::vector<std::string> &strformula,
                const std::vector<aminoacid::Aminoacids::Aminoacid> &aminoacids, const std::vector<int> &position,
                const std::string &graph, const std::string& extra_info);

        std::vector<Segment> containNRPsPart(nrpsprediction::NRPsPart predict_part) const override;

        NRPType getType() const override;
        std::vector<NRPLine> getLines() const override;

        bool is_valid_seg(int l, int r, int stp) const override;
    };
}


#endif //NRPSMATCHER_NRPLINE_H
