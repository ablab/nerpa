//
// Created by olga on 13.11.19.
//

#ifndef NERPA_ORDEREDGENESSCOREBASE_H
#define NERPA_ORDEREDGENESSCOREBASE_H

#include <Matcher/Score/Base/Score.h>

namespace matcher {
    class OrderedGenesScoreBase : public Score {
    private:
        double skip_segment_score = -1;
        double insertion = -1;
        double deletion = -1;
        double mismatch = -1;
        double open_gap = -1;
        double continue_gap = -1;

    public:
        explicit OrderedGenesScoreBase(std::unique_ptr<Score> base) : Score(std::move(base)) {}
        explicit OrderedGenesScoreBase(std::unique_ptr<Score> base,
                                       double skip_segment_score = -1,
                                       double insertion = -1,
                                       double deletion = -1,
                                       double mismatch = -1,
                                       double open_gap = -1,
                                       double continue_gap = -1) :
                Score(std::move(base)), skip_segment_score(skip_segment_score),
                insertion(insertion), deletion(deletion),
                mismatch(mismatch),
                open_gap(open_gap), continue_gap(continue_gap) {}

        double minScore(const int len) const override {
            return -len;
        }

        double maxScore(const int len) const override {
            return len;
        }

        double resultScore(double score, const int len) const override {
            return score/maxScore(len);
        }

        double SkipSegment() const override {
            return skip_segment_score;
        }

        double Mismatch() const override {
            return mismatch;
        }

        double openGap() const override {
            return open_gap;
        }

        double continueGap() const override {
            return continue_gap;
        }

        double InsertionScore() const override {
            return insertion;
        }

        double DeletionScore() const override {
            return deletion;
        }
    };
}


#endif //NERPA_ORDEREDGENESSCOREBASE_H
