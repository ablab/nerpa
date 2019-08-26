#ifndef NRPSMATCHER_NRPTAIL_H
#define NRPSMATCHER_NRPTAIL_H

#include <memory>
#include "NRP.h"
#include "NRPLine.h"


namespace nrp {
    class NRPtail : public NRP {
    private:
        std::shared_ptr<nrp::NRP> v1;
        std::shared_ptr<nrp::NRP> v2;
    public:
        NRPtail(std::shared_ptr<nrp::NRP> v1, std::shared_ptr<nrp::NRP> v2): v1(v1), v2(v2) {
        }

        NRPType getType() const override;

        int getLen() const override;

        int getInd(int i) const override;

        std::string getFormula(int i) const override;

        aminoacid::Aminoacid getAminoacid(int i) const override;

        void print() const override;

        std::string getGraphInString() const override;

        std::string get_file_name() const override;

        std::string get_extra_info() const override;

        std::vector<std::shared_ptr<NRP>> getLines() const override;

        bool is_valid_seg(int l, int r, int stp) const override;

        NRPtail(const NRP &refNrp);
    };
}


#endif //NRPSMATCHER_NRPTAIL_H
