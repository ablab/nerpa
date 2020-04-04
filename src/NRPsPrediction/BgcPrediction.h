#ifndef NRPSMATCHER_NRPSPREDICTION_H
#define NRPSMATCHER_NRPSPREDICTION_H

#include "OrfPrediction.h"
#include <iostream>

namespace nrpsprediction {
    //structure for store NRPs predictions
    class BgcPrediction {
    private:
        std::vector<OrfPrediction> orfs_;
        std::vector<OrfPrediction> short_orfs_;
    public:
        BgcPrediction() = default;
        BgcPrediction(std::vector<OrfPrediction> orfs): orfs_(orfs) {}
        BgcPrediction(std::vector<OrfPrediction> orfs, std::vector<OrfPrediction> short_orfs):
                orfs_(orfs), short_orfs_(short_orfs) {}
        const std::vector<OrfPrediction>& getOrfs() const;
        const std::vector<OrfPrediction>& getShortOrfs() const;

        int getSumPredictionLen() const;
    };
}


#endif //NRPSMATCHER_NRPSPREDICTION_H
