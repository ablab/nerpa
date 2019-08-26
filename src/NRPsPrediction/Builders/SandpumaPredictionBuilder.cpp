//
// Created by olga on 05.03.19.
//

#include <fstream>
#include <sstream>
#include <algorithm>
#include "SandpumaPredictionBuilder.h"

namespace nrpsprediction {
    void SandpumaPredictionBuilder::read_file(std::string file_name) {
        std::ifstream in(file_name);
        Token token;
        std::vector<Token> tokens;
        while (parse_token(in, token)) {
            tokens.push_back(token);
            token.res.resize(0);
        }

        std::sort(tokens.begin(), tokens.end());
        std::cerr << tokens.size() << "\n";
        int cur_id = 1;
        for (Token token : tokens) {
            if (!nrpparts.empty() && nrpparts.back().get_orf_name() == token.orf_id) {
                nrpparts[nrpparts.size() - 1].add_prediction(cur_id,
                                                             AminoacidPrediction(cur_id, parse_predictions(
                                                                     token)));
                ++cur_id;
            } else {
                cur_id = 1;
                if (nrpparts.size() > 0 && nrpparts[nrpparts.size() - 1].getAminoacidsPrediction().size() < 2) {
                    short_parts.push_back(nrpparts.back());
                    nrpparts.pop_back();
                }
                nrpparts.push_back(NRPsPart(file_name, token.orf_id, cur_id,
                                            AminoacidPrediction(cur_id,
                                                                parse_predictions(token))));
                ++cur_id;
            }
        }
        if (nrpparts.size() > 0 && nrpparts[nrpparts.size() - 1].getAminoacidsPrediction().size() < 2) {
            short_parts.push_back(nrpparts.back());
            nrpparts.pop_back();
        }
    }

    bool SandpumaPredictionBuilder::parse_token(std::ifstream &in, SandpumaPredictionBuilder::Token &token) {
        std::string s;
        if (getline(in, s)) {
            for (int i = 0; i < 5; ++i) {
                std::stringstream ss(s);
                std::string name;
                std::string source;
                std::string pred;
                ss >> name >> source >> pred;
                auto orf_id = split_name(name);
                token.orf_id = orf_id.first;
                token.id = orf_id.second;
                token.res.push_back(pred);
                if (i < 4) {
                    getline(in, s);
                }
            }
            return true;
        } else {
            return false;
        }
    }

    std::pair<std::string, int> SandpumaPredictionBuilder::split_name(std::string name) {
        std::stringstream ss(name);
        std::vector<std::string> split;
        std::string s;
        while (std::getline(ss, s, '_')) {
            split.push_back(s);
        }
        int id = 0;
        for (int i = 1; i < split[split.size() - 1].size(); ++i) {
            id = id * 10 + (split[split.size() - 1][i] - '0');
        }

        return std::make_pair(split[split.size() - 3] + "_" + split[split.size() - 2], id);
    }

    std::vector<AminoacidPrediction::AminoacidProb>
    SandpumaPredictionBuilder::parse_predictions(SandpumaPredictionBuilder::Token &t) {
        std::map<std::string, int> aa_score;
        for (std::string sres : t.res) {
            std::stringstream ss(sres);
            std::string aa;
            while (std::getline(ss, aa, '|')) {
                if (aa_score.count(aa) == 0) {
                    aa_score[aa] = 0;
                }
                aa_score[aa] += 30;
            }
        }

        std::vector<AminoacidPrediction::AminoacidProb> aminoacid_prediction;

        for (auto& item : aa_score) {
            if (getAAbyName(item.first) != aminoacid::AminoacidInfo::AMINOACID_CNT - 1) {
                aminoacid_prediction.push_back(AminoacidPrediction::AminoacidProb(
                        aminoacid::Aminoacid(getAAbyName(item.first)), item.second));
            }
        }

        std::sort(aminoacid_prediction.begin(), aminoacid_prediction.end(), [](AminoacidPrediction::AminoacidProb a,
                                                                               AminoacidPrediction::AminoacidProb b) {
            return a.prob > b.prob;
        });

        return aminoacid_prediction;
    }
}