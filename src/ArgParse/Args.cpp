//
// Created by olga on 19.07.19.
//

#include <fstream>
#include "Args.h"

Args::Args(std::string cfg_filename) {
    std::ifstream in(cfg_filename);
    in >> predictor_name;
    std::string tmp, val;
    in >> tmp >> val;
    insertion = (val == "on");
    in >> tmp >> val;
    deletion = (val == "on");
    double x;
    in >> tmp >> x;
    open_gap = x;
    in >> tmp >> x;
    continue_gap = x;
    in >> tmp >> val;
    single_match = (val == "on");
    in >> tmp >> x;
    single_match_coeff = x;
    modification = (val == "on");
    in >> modification_cfg;
    in >> AAmod_cfg;
}
