//
// Created by olga on 22.01.19.
//

#ifndef NRPSMATCHER_SCORE_H
#define NRPSMATCHER_SCORE_H

#include <NRP/NRP.h>
#include <Matcher/Segment.h>
#include <memory>
#include <map>
#include <algorithm>

namespace matcher {
    class Score {
    protected:
        double insertion = -1;
        double deletion = -5;

        std::map<double, double, std::greater<>> ProbGenCorrect;
        std::map<double, double, std::greater<>> ProbGenIncorrect;
//        //[100, 90, 80, 70, <= 60]
//        std::vector<double> ProbGenCorrect = {-0.07, -0.16, -0.7, -1.6, -1.1};
//        //[100, 90, 80, 70, <= 60]
//        std::vector<double> ProbGenIncorrect = {-2.66, -1.90, -0.7, -0.22, 0.};
//        //[100, 90, 80, 70, <=60]
//        std::vector<double> ProbGetScore = {-0.52, -2.06, -2.35, -2.14, -2.68};

        double probGenAA(const aminoacid::Aminoacid &nrpAA) const;
    public:
        explicit Score(const std::string &config);

        explicit Score(double insertion = -1, double deletion = -1) {
            this->insertion = insertion;
            this->deletion = deletion;
            ProbGenCorrect[0] = 1;
            ProbGenIncorrect[0] = -1;

            baseScore = nullptr;
            double curscore = 1;
            for (int i = 0; i < 100; ++i) {
                posscore[i] = curscore;
                curscore /= 1.25;
            }
        }

        explicit Score(std::unique_ptr<Score> base);

        virtual double minScore(const int len) const {
            if (baseScore != nullptr) {
                return baseScore->minScore(len);
            } else {
                return -len;
            }
        }

        virtual double resultScore(double score, const int len,
                                   const std::vector<Segment>& matched_parts,
                                   const nrpsprediction::BgcPrediction& prediction,
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

        virtual bool getScoreForSegment(const std::vector<aminoacid::Aminoacid>& amns,
                                        const nrpsprediction::BgcPrediction& prediction, int part_id, double& score) const;

        virtual double aaScore(const nrpsprediction::AAdomainPrediction &apred,
                       const aminoacid::Aminoacid &aminoacid) const;

        virtual std::pair<double, aminoacid::Aminoacid> getTheBestAAInPred(const nrpsprediction::AAdomainPrediction &apred,
                                                                           const aminoacid::Aminoacid &aminoacid,
                                                                           nrpsprediction::AAdomainPrediction::AminoacidProb &probRes,
                                                                           std::pair<int, int> &posRes) const;

        //return true if match is possible, false for mismatch
        virtual bool getScore(const aminoacid::Aminoacid& nrpAA,
                                const aminoacid::Aminoacid& predAA,
                                const nrpsprediction::AAdomainPrediction::AminoacidProb& prob,
                                const std::pair<int, int>& pos,
                                double& score) const;

        virtual double InsertionScore() const {
            if (baseScore != nullptr) {
                return baseScore->InsertionScore();
            } else {
                return insertion;
            }
        }

        virtual double DeletionScore() const {
            if (baseScore != nullptr) {
                return baseScore->DeletionScore();
            } else {
                return deletion;
            }
        }

        virtual double Mismatch(const aminoacid::Aminoacid& structure_aa, const nrpsprediction::AAdomainPrediction& aa_prediction) const;

    protected:
        std::unique_ptr<Score> baseScore;
        double posscore[100]{};
    };
}


#endif //NRPSMATCHER_SCORE_H
