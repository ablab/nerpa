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

        Match isCover(nrpsprediction::NRPsPrediction nrPsPrediction) override;

        NRPType getType() override;

        std::vector<Segment> containNRPsPart(nrpsprediction::NRPsPart predict_part) override;

        int getLen() override;

        int getInd(int i) override;

        std::string getFormula(int i) override;

        aminoacid::Aminoacids::Aminoacid getAminoacid(int i) override;

        void print() override;

        std::string getGraphInString() override;

        std::string get_file_name() override;

        std::string get_extra_info() override;
    };
}


#endif //NRPSMATCHER_NRPTAIL_H
