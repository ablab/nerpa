//
// Created by olga on 19.07.19.
//

#include <fstream>
#include "Args.h"

Args::Args(std::string &cfg_filename) {

    // TODO: make config independent of the key order (consider using some lightweight parser)
    std::ifstream in(cfg_filename);
    std::string tmp, val;
    double x;
    in >> tmp >> x;
    insertion = x;
    in >> tmp >> x;
    deletion = x;
    in >> modification_cfg;
    in >> monomer_cfg;
    in >> monomer_logP_cfg;
    in >> tmp >> threads;
    in >> tmp >> min_score;
    in >> tmp >> min_explain_part;
}
