//
// Created by olga on 23.02.19.
//

#include <fstream>
#include <sstream>
#include "MinowaPredictionBuilder.h"

namespace nrpsprediction {
    const std::string MinowaPredictionBuilder::AMINOACID_NAMES[aminoacid::Aminoacid::AMINOACID_CNT] = {"Trp", "Ser", "Gly", "uda", "Thr",
                                                                              "dhp", "Gln", "Dab", "Arg", "Lys",
                                                                              "ala-d", "Phe", "Val", "cha", "DHpg",
                                                                              "phg", "His", "Aeo", "Bmt", "hse",
                                                                              "met", "Ala", "tcl", "Sal", "allothr",
                                                                              "B-Ala", "dhb", "Ile", "end", "Leu",
                                                                              "gua", "homoTyr", "Glu", "bht", "Hpg",
                                                                              "apa", "Pro", "Tyr", "hyv", "Asn",
                                                                              "cit", "vol", "Cys", "Asp", "dht",
                                                                              "Ahp", "Orn", "apc", "Abu", "Aad",
                                                                              "pipecolate", "dpg", "none"};


    void MinowaPredictionBuilder::read_file(std::string file_name) {
        std::ifstream in(file_name);
        std::string s;
        while (getline(in, s)) {
            while (s != "\\\\" && getline(in, s)) {}
            std::string orf_name;
            getline(in, orf_name);
            getline(in, s);
            std::pair <std::string, int> orf_name_num = get_orf_name_and_order(orf_name);
            auto prediction = parse_predictions(in);
            if (!nrpparts.empty() && nrpparts[nrpparts.size() - 1].get_orf_name() == orf_name_num.first) {
                nrpparts[nrpparts.size() - 1].add_prediction(orf_name_num.second,
                                                             AminoacidPrediction(orf_name_num.second, prediction));
            } else {
                if (nrpparts.size() > 0 && nrpparts[nrpparts.size() - 1].getAminoacidsPrediction().size() < 2) {
                    nrpparts.pop_back();
                }
                nrpparts.push_back(NRPsPart(file_name, orf_name_num.first, orf_name_num.second,
                                            AminoacidPrediction(orf_name_num.second, prediction)));
            }
        }

        if (nrpparts.size() > 0 && nrpparts[nrpparts.size() - 1].getAminoacidsPrediction().size() < 2) {
            nrpparts.pop_back();
        }
        /*std::cerr << nrpparts.size() << "\n";
        for (int i = 0; i < nrpparts.size(); ++i) {
            auto predictions = nrpparts[i].getAminoacidsPrediction();
            std::cerr << i << " size:" << predictions.size() << "\n";
            for (int j = 0; j < predictions.size(); ++j) {
                auto probs = predictions[j].getAAPrediction();
                for (int g = 0; g < probs.size(); ++g) {
                    std::cerr << probs[g].aminoacid << " " << probs[g].prob << "; ";
                }
                std::cerr << "\n";
            }
        }*/
        in.close();
    }

    std::vector<AminoacidPrediction::AminoacidProb>
    MinowaPredictionBuilder::parse_predictions(std::ifstream &in) {
        std::string name;
        double score;

        std::vector<std::pair<std::string, double> > aacids;
        std::string s;
        getline(in, s);
        while (s != "") {
            std::stringstream ss(s);
            ss >> name >> score;
            aacids.push_back(std::make_pair(name, score));
            getline(in, s);
        }

        std::vector<AminoacidPrediction::AminoacidProb> aminoacid_prediction;


        double val = std::min(aacids[0].second, std::max(aacids[2].second, 100.));

        for (int i = 0; i < aacids.size(); ++i) {
            if (aacids[i].second >= val - EPS) {
                aminoacid_prediction.push_back(AminoacidPrediction::AminoacidProb(
                        aminoacid::Aminoacid(getAAbyName(aacids[i].first, AMINOACID_NAMES)), aacids[i].second));
            }
        }
        return aminoacid_prediction;
    }
}
