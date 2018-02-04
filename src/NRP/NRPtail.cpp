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

std::vector<nrp::NRP::Segment> nrp::NRPtail::containNRPsPart(nrpsprediction::NRPsPart predict_part) {
    std::vector<nrp::NRP::Segment> seg1 = v1.containNRPsPart(predict_part), seg2 = v2.containNRPsPart(predict_part);

    seg1.insert(seg1.end(), seg2.begin(), seg2.end());
    return seg1;
}
