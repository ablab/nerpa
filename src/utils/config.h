#pragma once

#include "cxxopts.hpp"

struct nerpa_config {
    std::string modification_cfg;
    std::string monomer_cfg;
    std::string monomer_logP_cfg;
    std::string prob_cfg;

    unsigned int threads;
    double min_score;
    double min_explain_part;
    double monomer_info_default_logp;
    double aminoacid_info_default_logp;

    bool debug;
};

void load(nerpa_config &cfg, const cxxopts::Options& options);

// config singleton-wrap
struct config {
    static void create_instance(const cxxopts::Options& options) {
        load(inner_cfg(), options);
        is_initialized() = true;
    }

    static nerpa_config const &get() {
        //VERIFY_MSG(is_initialized(), "Config not initialized");
        return inner_cfg();
    }

    static nerpa_config &get_writable() {
        //VERIFY_MSG(is_initialized(), "Config not initialized");
        return inner_cfg();
    }

private:
    static nerpa_config &inner_cfg() {
        static nerpa_config config;
        return config;
    }

    static bool &is_initialized() {
        static bool is_initialized = false;
        return is_initialized;
    }
};
