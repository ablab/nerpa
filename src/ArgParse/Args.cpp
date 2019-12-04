//
// Created by olga on 19.07.19.
//

#include <fstream>
#include "Args.h"

Args::Args(std::string cfg_filename) {
    std::ifstream in(cfg_filename);
    in >> predictor_name;
    std::string tmp, val;
    double x;
    in >> tmp >> x;
    insertion = x;
    in >> tmp >> x;
    deletion = x;
    in >> tmp >> x;
    open_gap = x;
    in >> tmp >> x;
    continue_gap = x;
    in >> tmp >> x;
    mismatch = x;
    in >> tmp >> x;
    skip_segment = x;
    in >> tmp >> val;
    modification = (val == "on");
    in >> modification_cfg;
    in >> AAmod_cfg;
}
