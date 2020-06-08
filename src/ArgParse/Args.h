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
    double open_gap;
    double continue_gap;
    double mismatch;
    double skip_segment;
    bool modification;
    std::string modification_cfg;
    std::string AAmod_cfg;
    std::string monomer_cfg;
    std::string monomer_logP_cfg;
    unsigned int threads;

    Args(std::string cfg_filename);
};


#endif //NRPSMATCHER_ARGS_H
