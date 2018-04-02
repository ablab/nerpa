#include <sstream>
#include <iostream>
#include "NRPsPredictionGenerator.h"

using namespace nrpsprediction;
typedef aminoacid::Aminoacids::Aminoacid AA;

namespace nrp_generator {
    NRPsPredictionGenerator::NRPsPredictionGenerator(
            std::vector<NRPsPrediction> &predictions) {
        for (int i = 0; i < 11; ++i) {
            cntAA[i].resize(aminoacid::Aminoacids::AMINOACID_CNT);
        }

        for (NRPsPrediction &pred : predictions) {
            std::vector<NRPsPart> parts = pred.getNrpsParts();
            for (NRPsPart part : parts) {
                std::vector<AminoacidPrediction> pred = part.getAminoacidsPrediction();
                for (int i = 0; i < pred.size(); ++i) {
                    std::vector<AminoacidPrediction::AminoacidProb> probs = pred[i].getAAPrediction();
                    for (auto &prob : probs) {
                        cntAA[int(prob.prob / 10)][int(prob.aminoacid)] += 1;
                    }
                }
            }
        }

        for (int i = 0; i < 11; ++i) {
            for (int j = 1; j < cntAA[i].size(); ++j) {
                cntAA[i][j] += cntAA[i][j - 1];
            }
        }
    }

    NRPsPrediction NRPsPredictionGenerator::genPrediction(NRPsPrediction prediction) {
        std::vector<NRPsPart> genparts;
        std::vector<NRPsPart> parts = prediction.getNrpsParts();

        for (NRPsPart part : parts) {
            NRPsPart newPart("", "");
            std::vector<AminoacidPrediction> pred = part.getAminoacidsPrediction();
            for (int i = 0; i < pred.size(); ++i) {
                std::stringstream ss;
                std::vector<AminoacidPrediction::AminoacidProb> probs = pred[i].getAAPrediction();
                std::vector<AA> gena;
                for (auto &prob : probs) {
                    ss << aminoacid::Aminoacids::AMINOACID_NAMES[genAA(prob.prob, gena)] << "(" << prob.prob << ");";
                }
                newPart.add_prediction(i + 1, ss.str());
            }
            genparts.push_back(newPart);
        }

        return nrpsprediction::NRPsPrediction(genparts);
    }

    AA NRPsPredictionGenerator::genAA(double prob, std::vector<AA> &aa) {
        int pos = int(prob/10);
        AA cur = AA::end;
        bool cont = true;
        while (cont) {
            int rnd = rand()%(cntAA[pos][cntAA[pos].size() - 1] + 1);
            cur = AA(std::lower_bound(cntAA[pos].begin(), cntAA[pos].end(), rnd) - cntAA[pos].begin());
            cont = false;
            for (int i = 0; i < aa.size(); ++i) {
                if (aa[i] == cur) {
                    cont = true;
                }
            }
        }
        aa.push_back(cur);
        return cur;
    }
}