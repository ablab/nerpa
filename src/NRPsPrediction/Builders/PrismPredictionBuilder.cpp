//
// Created by olga on 01.03.19.
//

#include "PrismPredictionBuilder.h"
#include <fstream>


namespace nrpsprediction {
    const std::string PrismPredictionBuilder::AMINOACID_NAMES[aminoacid::Aminoacid::AMINOACID_CNT] = {"W", "S", "G", "uda", "T",
                                                                                                       "dhp", "Q", "Dab", "R", "K",
                                                                                                       "ala-d", "F", "V", "cha", "DHpg",
                                                                                                       "phg", "H", "Aeo", "Bmt", "hse",
                                                                                                       "M", "A", "tcl", "Sal", "allothr",
                                                                                                       "B-Ala", "Dhab", "I", "end", "L",
                                                                                                       "gua", "homoTyr", "E", "bht", "Hpg",
                                                                                                       "apa", "P", "Y", "hyv", "N",
                                                                                                       "cit", "vol", "C", "D", "dht",
                                                                                                       "Ahp", "Orn", "apc", "Abu", "Aad",
                                                                                                       "Pip", "dpg", "none"};

    void PrismPredictionBuilder::read_file(std::string file_name) {
        using json = nlohmann::json;

        json jsn;
        std::ifstream in(file_name);
        in >> jsn;

        for (auto cluster : jsn["prism_results"]["clusters"]) {
            for (auto orf : cluster["orfs"]) {
                if (orf["type"] != "nrps") {
                    continue;
                }

                std::string name = orf["name"];
                int curId = 0;
                for (int i = 0; i < orf["domains"].size(); ++i) {
                    if (orf["domains"][i]["name"] == "A") {
                        if (curId > 0) {
                            nrpparts.back().add_prediction(curId + 1,
                                                           AminoacidPrediction(curId + 1,
                                                                               parse_predictions(
                                                                                       orf["domains"][i]["substrates"])));
                        } else {
                            if (nrpparts.size() > 0 && nrpparts.back().getAminoacidsPrediction().size() < 2) {
                                nrpparts.pop_back();
                            }
                            nrpparts.push_back(NRPsPart(file_name, name, curId + 1,
                                                        AminoacidPrediction(curId + 1,
                                                                            parse_predictions(
                                                                                    orf["domains"][i]["substrates"]))));

                        }

                        curId += 1;
                    }
                }
            }
        }

        if (nrpparts.size() > 0 && nrpparts.back().getAminoacidsPrediction().size() < 2) {
            nrpparts.pop_back();
        }
        in.close();
    }

    std::vector<AminoacidPrediction::AminoacidProb>
    PrismPredictionBuilder::parse_predictions(nlohmann::json predictions) {
        std::vector<AminoacidPrediction::AminoacidProb> aminoacid_prediction;

        for (auto aa : predictions) {
            aminoacid_prediction.push_back(AminoacidPrediction::AminoacidProb(
                    aminoacid::Aminoacid(getAAbyName(aa["name"], AMINOACID_NAMES)), aa["score"]));
        }

        return aminoacid_prediction;
    }
}