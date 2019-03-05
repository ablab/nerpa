//
// Created by olga on 05.03.19.
//

#include <fstream>
#include <sstream>
#include <algorithm>
#include "SandpumaPredictionBuilder.h"

namespace nrpsprediction {
    const std::string SandpumaPredictionBuilder::AMINOACID_NAMES[aminoacid::Aminoacid::AMINOACID_CNT] = {"trp", "ser", "gly", "uda", "thr",
                                                                              "dhp", "gln", "dab", "arg", "lys",
                                                                              "ala-d", "phe", "val", "cha", "dhpg",
                                                                              "phg", "his", "aeo", "bmt", "hse",
                                                                              "met", "ala", "tcl", "sal", "allothr",
                                                                              "b-ala", "dhb", "ile", "end", "leu",
                                                                              "gua", "hty", "glu", "bht", "hpg",
                                                                              "apa", "pro", "tyr", "hyv", "asn",
                                                                              "cit", "vol", "cys", "asp", "dht",
                                                                              "ahp", "orn", "apc", "abu", "aad",
                                                                              "pip", "dpg", "none"};

    void SandpumaPredictionBuilder::read_file(std::string file_name) {
        file_name += "/ind.res.tsv";
        std::ifstream in(file_name);
        Token token;
        std::vector<Token> tokens;
        while (parse_token(in, token)) {
            tokens.push_back(token);
        }

        std::sort(tokens.begin(), tokens.end());
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
                    nrpparts.pop_back();
                }
                nrpparts.push_back(NRPsPart(file_name, token.orf_id, cur_id,
                                            AminoacidPrediction(cur_id,
                                                                parse_predictions(token))));
                ++cur_id;
            }
        }
        if (nrpparts.size() > 0 && nrpparts[nrpparts.size() - 1].getAminoacidsPrediction().size() < 2) {
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

        return std::make_pair(split[split.size() - 3], std::stoi(split.back()));
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
                aa_score[aa] += 20;
            }
        }

        std::vector<AminoacidPrediction::AminoacidProb> aminoacid_prediction;

        for (auto& item : aa_score) {
            aminoacid_prediction.push_back(AminoacidPrediction::AminoacidProb(
                    aminoacid::Aminoacid(getAAbyName(item.first, AMINOACID_NAMES)), item.second));
        }

        return std::vector<AminoacidPrediction::AminoacidProb>();
    }
}