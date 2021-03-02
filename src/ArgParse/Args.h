//
// Created by olga on 19.07.19.
//

#ifndef NRPSMATCHER_ARGS_H
#define NRPSMATCHER_ARGS_H

#include <string>

class Args {
public:
    double insertion;
    double deletion;
    std::string modification_cfg;
    std::string monomer_cfg;
    std::string monomer_logP_cfg;
    unsigned int threads;
    double min_score;
    double min_explain_part;

    explicit Args(std::string &cfg_filename);
};


#endif //NRPSMATCHER_ARGS_H
