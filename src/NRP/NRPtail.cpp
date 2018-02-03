#include <iostream>
#include "NRPtail.h"

nrp::NRP::Match nrp::NRPtail::isCover(nrpsprediction::NRPsPrediction nrPsPrediction) {
    Match m1 = v1.isCover(nrPsPrediction);
    Match m2 = v2.isCover(nrPsPrediction);
    if (m1.score() > m2.score()) {
        return m1;
    } else {
        return m2;
    }
}
