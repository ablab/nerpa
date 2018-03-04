#ifndef NRPSMATCHER_NORMALIZEDMATCH_H
#define NRPSMATCHER_NORMALIZEDMATCH_H

#include <NRP/NRP.h>
#include <NRP/NRPGenerator.h>

namespace normalized_match {
    using namespace nrp;
    class NormalizedMatch {
    private:
        const int CNT_GEN = 100;
        double score = 0;

        NRP::Match match;
    public:
        NormalizedMatch(NRP::Match match, NRPGenerator generator, nrpsprediction::NRPsPrediction prediction, NRP* mol);
        void print(std::ofstream& out);
        void print_short(std::ofstream& out);
        void print_short_prediction(std::ofstream& out);

        double calcMean(std::vector<int> score);

        double calcSD(std::vector<int> score, double mean);

        bool operator < (NormalizedMatch b);
    };
}


#endif //NRPSMATCHER_NORMALIZEDMATCH_H
