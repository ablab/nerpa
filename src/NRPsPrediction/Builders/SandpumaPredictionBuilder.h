//
// Created by olga on 05.03.19.
//

#ifndef NRPSMATCHER_SANDPUMAPREDICTIONBUILDER_H
#define NRPSMATCHER_SANDPUMAPREDICTIONBUILDER_H

#include "PredictionBuilderBase.h"

namespace nrpsprediction {
    class SandpumaPredictionBuilder : public PredictionBuilderBase {
    private:
        struct Token {
            std::string orf_id;
            int id;
            std::vector<std::string> res;

            bool operator < (const Token &t) const {
                return orf_id < t.orf_id || (orf_id == t.orf_id && id < t.id);
            }
        };

        bool parse_token(std::ifstream& in, Token& t);

        std::vector<AAdomain_Prediction::AminoacidProb> parse_predictions(Token& t);
    public:
        //parse ctg1_minowa_nrpspredoutput.txt file
        //expected groupded by orfs and sorted by num in one group
        void read_file(std::string file_name) override;

        std::pair<std::string, int> split_name(std::string name);
    };
}


#endif //NRPSMATCHER_SANDPUMAPREDICTIONBUILDER_H
