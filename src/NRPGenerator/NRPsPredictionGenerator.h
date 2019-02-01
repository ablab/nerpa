#ifndef NRPSMATCHER_NRPSPREDICTIONGENERATOR_H
#define NRPSMATCHER_NRPSPREDICTIONGENERATOR_H


#include <vector>
#include <NRPsPrediction/NRPsPrediction.h>

namespace nrp_generator {
    class NRPsPredictionGenerator {
    private:
        std::vector<int> cntAA[11];
        aminoacid::Aminoacid::AminoacidId genAA(double prob, std::vector<aminoacid::Aminoacid::AminoacidId>& aa);
    public:
        explicit NRPsPredictionGenerator(std::vector<nrpsprediction::NRPsPrediction>& predictions);
        nrpsprediction::NRPsPrediction genPrediction(nrpsprediction::NRPsPrediction prediction);
    };
}

#endif //NRPSMATCHER_NRPSPREDICTIONGENERATOR_H
