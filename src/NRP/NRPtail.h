#ifndef NRPSMATCHER_NRPTAIL_H
#define NRPSMATCHER_NRPTAIL_H

#include "NRP.h"
#include "NRPLine.h"


namespace nrp {
    class NRPtail : public NRP {
    private:
        NRPLine v1;
        NRPLine v2;
    public:
        NRPtail(NRPLine v1, NRPLine v2): v1(v1), v2(v2) {
        }

        NRPType getType() const override;

        std::vector<Segment> containNRPsPart(nrpsprediction::NRPsPart predict_part) const override;

        int getLen() const override;

        int getInd(int i) const override;

        std::string getFormula(int i) const override;

        aminoacid::Aminoacids::Aminoacid getAminoacid(int i) const override;

        void print() const override;

        std::string getGraphInString() const override;

        std::string get_file_name() const override;

        std::string get_extra_info() const override;

        std::vector<NRPLine> getLines() const override;

        bool is_valid_seg(int l, int r, int stp) const override;
    };
}


#endif //NRPSMATCHER_NRPTAIL_H
