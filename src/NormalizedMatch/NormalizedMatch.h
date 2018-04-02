#ifndef NRPSMATCHER_NORMALIZEDMATCH_H
#define NRPSMATCHER_NORMALIZEDMATCH_H

#include <NRP/NRP.h>
#include <NRPGenerator/NRPGenerator.h>
#include <NRPGenerator/NRPsPredictionGenerator.h>

namespace normalized_match {
    using namespace nrp;
    class NormalizedMatch {
    private:
        static const int CNT_GEN;
        static const double EPS;
        double score = 0;
        double p_value = 0;

        NRP::Match match;
    public:
        NormalizedMatch() = default;
        NormalizedMatch(NRP::Match match, nrp_generator::NRPGenerator* generator,
                        nrpsprediction::NRPsPrediction prediction, NRP* mol);

        NormalizedMatch(NRP::Match match, nrp_generator::NRPsPredictionGenerator* generator,
                        nrpsprediction::NRPsPrediction prediction, NRP* mol);
        void print(std::ofstream& out);
        void print_short(std::ofstream& out);
        void print_short_prediction(std::ofstream& out);
        void print_csv(std::ofstream& out);

        double calcMean(std::vector<double> score);

        double calcSD(std::vector<double> score, double mean);

        bool operator < (const NormalizedMatch &b) const;
    };
}


#endif //NRPSMATCHER_NORMALIZEDMATCH_H
