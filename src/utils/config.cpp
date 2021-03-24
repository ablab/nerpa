#include <algorithm>
#include "config.h"

void load(nerpa_config &cfg, const cxxopts::Options& options) {
    std::string configs_dir = options["configs_dir"].as<std::string>() + "/";

    cfg.modification_cfg = configs_dir + "modifications.tsv";
    cfg.monomer_cfg = configs_dir + "monomers.tsv";
    cfg.monomer_logP_cfg = configs_dir + "monomersLogP.tsv";
    cfg.prob_cfg = configs_dir + "prob_gen.cfg";

    cfg.threads = std::max((unsigned int)1, options["threads"].as<unsigned int>());
    cfg.min_score = options["min_score"].as<double>();
    cfg.min_explain_part = options["min_explain_part"].as<double>();
    cfg.monomer_info_default_logp = options["default_monomer_logp"].as<double>();
    cfg.aminoacid_info_default_logp = options["default_aminoacid_logp"].as<double>();

    cfg.debug = false;
}