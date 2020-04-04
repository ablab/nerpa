//
// Created by olga on 01.03.19.
//

#include "PrismPredictionBuilder.h"
#include <fstream>


namespace nrpsprediction {
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
                                                           AAdomainPrediction(curId + 1,
                                                                              parse_predictions(
                                                                                       orf["domains"][i]["substrates"])));
                        } else {
                            if (nrpparts.size() > 0 && nrpparts.back().getAAdomainPrediction().size() < 2) {
                                nrpparts.pop_back();
                            }
                            nrpparts.push_back(OrfPrediction(file_name, name, curId + 1,
                                                             AAdomainPrediction(curId + 1,
                                                                                 parse_predictions(
                                                                                    orf["domains"][i]["substrates"]))));

                        }

                        curId += 1;
                    }
                }
            }
        }

        if (nrpparts.size() > 0 && nrpparts.back().getAAdomainPrediction().size() < 2) {
            nrpparts.pop_back();
        }
        in.close();
    }

    std::vector<AAdomainPrediction::AminoacidProb>
    PrismPredictionBuilder::parse_predictions(nlohmann::json predictions) {
        std::vector<AAdomainPrediction::AminoacidProb> aminoacid_prediction;

        for (auto aa : predictions) {
            aminoacid_prediction.push_back(AAdomainPrediction::AminoacidProb(
                    aminoacid::Aminoacid(getAAbyName(aa["name"])), aa["score"]));
        }

        return aminoacid_prediction;
    }
}