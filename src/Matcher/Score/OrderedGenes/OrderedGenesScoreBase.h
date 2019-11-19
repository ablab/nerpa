//
// Created by olga on 13.11.19.
//

#ifndef NERPA_ORDEREDGENESSCOREBASE_H
#define NERPA_ORDEREDGENESSCOREBASE_H

#include <Matcher/Score/Base/Score.h>

namespace matcher {
    class OrderedGenesScoreBase : public Score {
    public:
        explicit OrderedGenesScoreBase(std::unique_ptr<Score> base) : Score(std::move(base)) {}

        double minScore(const int len) const override {
            return -len;
        }

        double maxScore(const int len) const override {
            return len;
        }
    };
}


#endif //NERPA_ORDEREDGENESSCOREBASE_H
