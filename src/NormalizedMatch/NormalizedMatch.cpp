#include <cmath>
#include "NormalizedMatch.h"

normalized_match::NormalizedMatch::NormalizedMatch(nrp::NRP::Match match, nrp::NRPGenerator generator,
                                                   nrpsprediction::NRPsPrediction prediction, nrp::NRP *mol) {
    this->match = match;
    std::vector<int> scores;
    for (int i = 0; i < CNT_GEN; ++i) {
        nrp::NRP* genNRP = generator.generate(mol->getType(), mol->getLen());
        scores.push_back(genNRP->isCover(prediction).score());
        delete genNRP;
    }

    double mean = calcMean(scores);
    double SD = calcSD(scores, mean);

    score = (match.score() - mean)/std::sqrt(SD);
}

double normalized_match::NormalizedMatch::calcMean(std::vector<int> score) {
    int sum = 0;
    for (int i = 0; i < score.size(); ++i) {
        sum += score[i];
    }

    return double(sum)/score.size();
}

double normalized_match::NormalizedMatch::calcSD(std::vector<int> score, double mean) {
    double sum = 0;
    for (int i = 0; i < score.size(); ++i) {
        sum += (score[i] - mean) * (score[i] - mean);
    }
    return sum/score.size();
}

bool normalized_match::NormalizedMatch::operator<(normalized_match::NormalizedMatch b) {
    return score > b.score;
}

void normalized_match::NormalizedMatch::print(std::ofstream &out) {
    match.print(out, score);
}

void normalized_match::NormalizedMatch::print_short(std::ofstream &out) {
    match.print_short(out, score);
}

void normalized_match::NormalizedMatch::print_short_prediction(std::ofstream &out) {
    match.print_short_prediction(out, score);
}
