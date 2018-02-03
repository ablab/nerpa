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
    };
}


#endif //NRPSMATCHER_NRPTAIL_H
