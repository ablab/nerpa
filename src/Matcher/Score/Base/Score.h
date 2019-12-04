//
// Created by olga on 22.01.19.
//

#ifndef NRPSMATCHER_SCORE_H
#define NRPSMATCHER_SCORE_H

#include <NRP/NRP.h>
#include <Matcher/Segment.h>
#include <memory>

namespace matcher {
    class Score {
    public:
        Score();
        explicit Score(std::unique_ptr<Score> base);
        virtual double minScore(const int len) const {
            if (baseScore != nullptr) {
                return baseScore->minScore(len);
            } else {
                return -len - 1;
            }
        }

        virtual double maxScore(const int len) const  {
            if (baseScore != nullptr) {
                return baseScore->maxScore(len);
            } else {
                int MAX_PRED_LEN = 4;
                return len - (len + 1) / MAX_PRED_LEN;
            }
        }

        virtual double resultScore(double score, const int len,
                                   const std::vector<Segment>& matched_parts,
                                   const nrpsprediction::NRPsPrediction& prediction,
                                   const nrp::NRP& nrp) const {
            if (baseScore != nullptr) {
                return baseScore->resultScore(score, len, matched_parts, prediction, nrp);
            } else {
                return score;
            }
        }

        virtual double resultScore(double score, const int len) const {
            if (baseScore != nullptr) {
                return baseScore->resultScore(score, len);
            } else {
                return score;
            }
        }

        virtual double openGap() const {
            if (baseScore != nullptr) {
                return  baseScore->openGap();
            } else {
                return -1;
            }
        }

        virtual double continueGap() const {
            if (baseScore != nullptr) {
                return baseScore->continueGap();
            } else {
                return 0;
            }
        }

        virtual double addSegment(Segment seg) const {
            if (baseScore != nullptr) {
                return baseScore->addSegment(seg);
            } else {
                return seg.scor - 1;
            }
        }

        virtual bool getScoreForSegment(const std::vector<aminoacid::Aminoacid>& amns,
                                const nrpsprediction::NRPsPrediction& prediction, int part_id, double& score) const;

        virtual double aaScore(const nrpsprediction::AminoacidPrediction &apred,
                       const aminoacid::Aminoacid &aminoacid) const;

        virtual std::pair<double, aminoacid::Aminoacid> getTheBestAAInPred(const nrpsprediction::AminoacidPrediction &apred,
                                                                           const aminoacid::Aminoacid &aminoacid,
                                                                           nrpsprediction::AminoacidPrediction::AminoacidProb &probRes,
                                                                           std::pair<int, int> &posRes) const;

        virtual double singleUnitScore(const nrpsprediction::AminoacidPrediction &apred,
                                       const aminoacid::Aminoacid &aminoacid) const;

        virtual double getScore(const aminoacid::Aminoacid& nrpAA,
                                const aminoacid::Aminoacid& predAA,
                                const nrpsprediction::AminoacidPrediction::AminoacidProb& prob,
                                const std::pair<int, int>& pos) const;

        virtual double InDelScore(double score, const int len) const {
            return score - 1./len;
        }


        virtual double InsertionScore() const {
            if (baseScore != nullptr) {
                return baseScore->InsertionScore();
            } else {
                return -1;
            }
        }

        virtual double DeletionScore() const {
            if (baseScore != nullptr) {
                return baseScore->DeletionScore();
            } else {
                return -1;
            }
        }

        virtual double SkipSegment() const {
            if (baseScore != nullptr) {
                return baseScore->SkipSegment();
            } else {
                return -1;
            }
        }

        virtual double Mismatch() const {
            if (baseScore != nullptr) {
                return baseScore->Mismatch();
            } else {
                return -1;
            }
        }

    protected:
        std::unique_ptr<Score> baseScore;
        double posscore[100];

    };
}


#endif //NRPSMATCHER_SCORE_H
