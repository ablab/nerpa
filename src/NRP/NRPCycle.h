#ifndef NRPSMATCHER_NRPCYCLE_H
#define NRPSMATCHER_NRPCYCLE_H

#include "NRP.h"
#include "NRPLine.h"

namespace nrp {
    class NRPCycle : public NRP {
    public:
        NRPCycle(const std::string &file_name, const std::vector<std::string> &strformula,
                 const std::vector<aminoacid::Aminoacid> &aminoacids, const std::vector<int> &position,
                 const std::string &graph, const std::string& extra_info);

        NRPType getType() const override;

        std::vector<std::shared_ptr<NRP>> getLines() const override;

        bool is_valid_seg(int l, int r, int stp) const override;

        std::string structure_to_string() const override;

        explicit NRPCycle(const NRP &refNrp);
    };
}


#endif //NRPSMATCHER_NRPCYCLE_H
