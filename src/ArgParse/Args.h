//
// Created by olga on 19.07.19.
//

#ifndef NRPSMATCHER_ARGS_H
#define NRPSMATCHER_ARGS_H

#include <string>

class Args {
public:
    std::string predictor_name;
    bool insertion;
    bool deletion;
    double open_gap;
    double continue_gap;
    bool single_match;
    double single_match_coeff;
    bool modification;
    std::string modification_cfg;
    std::string AAmod_cfg;

    Args(std::string cfg_filename);
};


#endif //NRPSMATCHER_ARGS_H
