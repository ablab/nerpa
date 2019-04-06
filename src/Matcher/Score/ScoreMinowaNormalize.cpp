//
// Created by olga on 06.04.19.
//

#include "ScoreMinowaNormalize.h"

namespace matcher {
    double ScoreMinowaNormalize::resultScore(double score, const int len) const {
        return score / maxScore(len);
    }
};
