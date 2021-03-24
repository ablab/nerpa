#include <iostream>
#include <fstream>
#include <NRP/NRPBuilder.h>
#include <NRP/MonomericNRPBuilder.h>
#include <algorithm>
#include <cstring>
#include "NRP/NRP.h"
#include "NRPsPrediction/BgcPrediction.h"
#include <Logger/log_writers.hpp>
#include <NRPsPrediction/Builders/Nrpspredictor2Builder.h>
#include <ArgParse/Args.h>
#include <Aminoacid/ModificationInfo.h>
#include "utils/openmp_wrapper.h"
#include "utils/cxxopts.hpp"
#include <Matcher/OrderedGenesMatcher.h>
#include <Aminoacid/MonomerInfo.h>

cxxopts::Options parse_options(int argc, char **argv) {
    cxxopts::Options options(argv[0],
                             " <file_with_predictions> <file_with_structures> <file_with_aa> <config_file> - matching BGC predictions with NRPs");
    options.add_options()
            ("h,help", "Print help")
            ("p,predictions", "File listing paths to antiSMASH predictions", cxxopts::value<std::string>(), "FILE")
            ("s,structures", "File with parsed NRP structures", cxxopts::value<std::string>(), "FILE")
            ("a,amino_acids", "File with amino acids", cxxopts::value<std::string>(), "FILE")
            ("c,config", "Configuration file", cxxopts::value<std::string>(), "FILE");

    options.add_options("Advanced")
            ("start_from", "Starting molecule index in the list of NRP structures", cxxopts::value<int>()->default_value("0"), "N");

    const std::vector<std::string> all_groups({"", "Advanced"});

    try {
        options.parse_positional(std::vector<std::string>({"predictions", "structures", "amino_acids", "config"}));
        options.parse(argc, argv);
    } catch (const cxxopts::OptionException &e) {
        std::cout << "error parsing options: " << e.what() << std::endl;
        std::cout << options.help(all_groups) << std::endl;
        exit(1);
    }

    if (options.count("help")) {
        std::cout << options.help(all_groups) << std::endl;
        exit(0);
    }
    for (const auto &required : std::vector<std::pair<std::string, std::string>>{
            {"predictions",   "Predictions file"},
            {"structures",   "Structures file"},
            {"amino_acids", "AA file"},
            {"config",   "Config file"}}) {
        if (options.count(required.first) != 1) {
            std::cout << required.second << " should be specified" << std::endl << std::endl;
            std::cout << options.help(all_groups) << std::endl;
            exit(0);
        }
    }

    return options;
}

void getPredictor(nrpsprediction::PredictionBuilderBase*& predictionBuilder) {
    predictionBuilder = new nrpsprediction::Nrpspredictor2Builder();
}

std::string get_file_name(std::string cur_line) {
    std::string res = "";
    for (int i = 0; i < cur_line.size() && cur_line[i] != '\t'; ++i) {
        res += cur_line[i];
    }
    return res;
}

std::vector<nrpsprediction::BgcPrediction>  save_predictions(const std::string& file_name) {
    std::vector<nrpsprediction::BgcPrediction> preds;
    std::ifstream in_predictions_files(file_name);

    INFO(file_name);
    std::string cur_prediction_file;
    std::string cur_line;

    while(getline(in_predictions_files, cur_line)) {
        //std::stringstream ss(cur_line);
        //ss >> cur_prediction_file;
        //if (cur_line[0] != '/') {
        //    cur_line = getDir(file_name) + "/" + cur_line;
        //}
        INFO(cur_line);
        std::string info_file_name = get_file_name(cur_line);
        INFO(info_file_name)

        nrpsprediction::PredictionBuilderBase* nrPsPredictionBuilder;
        getPredictor(nrPsPredictionBuilder);
        nrPsPredictionBuilder->read_file(info_file_name);

        preds.push_back(nrPsPredictionBuilder->getPrediction());
        INFO("Parts in prediction: " << preds.back().getOrfs().size());
        delete(nrPsPredictionBuilder);
    }

    return preds;
}

std::vector<std::shared_ptr<nrp::NRP>> load_nrps_from_monomeric_info(const std::string& file_name) {
    std::vector<std::shared_ptr<nrp::NRP>> nrps;
    std::ifstream in_nrps_files(file_name);
    std::string cur_line;
    std::string cur_id;

    while(getline(in_nrps_files, cur_line)) {
        INFO(cur_line);
        std::stringstream ss(cur_line);
        std::string extra;
        ss >> cur_id;
        getline(ss, extra);
        try {
            std::shared_ptr<nrp::NRP> nrp_from_fragment_graph = nrp::MonomericNRPBuilder::build(cur_id, extra);
            if (nrp_from_fragment_graph == nullptr) {
                continue;
            }
            nrp_from_fragment_graph->print();
            nrps.push_back(nrp_from_fragment_graph);
        } catch (...) {
            continue;
        }
    }
    return nrps;
}

void getScoreFunction(Args args, matcher::Score*& score) {
    using namespace matcher;
    score = new Score(args.prob_cfg);
}

matcher::MatcherBase* getMatcher(Args args) {
    return new matcher::OrderedGenesMatcher();
}


void run_mol_predictions(std::vector<nrpsprediction::BgcPrediction> preds, std::shared_ptr<nrp::NRP> mol, std::string output_filename,
                         Args args) {
    matcher::Score* score;
    getScoreFunction(args, score);
    std::vector<matcher::MatcherBase::Match> nrpsMatchs;
    for (int i = 0; i < preds.size(); ++i) {
        if (preds[i].getOrfs().size() == 0) continue;
        //std::cerr << mol->get_file_name() << " " << preds[i].getOrfs()[0].get_file_name() << "\n";
        matcher::MatcherBase* matcher = getMatcher(args);
        matcher::MatcherBase::Match match = matcher->getMatch(mol, &preds[i], score);
        delete matcher;

        if (match.score() >= args.min_score &&
                (double)match.getCntMatch()/preds[i].getSumPredictionLen() >= args.min_explain_part) {
            nrpsMatchs.push_back(match);
            std::ofstream out(output_filename);
            match.print(out);
        }

    }

    INFO("Found: " << nrpsMatchs.size() << " predictions");

    std::sort(nrpsMatchs.begin(), nrpsMatchs.end());
    if (nrpsMatchs.size() == 0) {
        return;
    }

#pragma omp critical
    {
        std::ofstream out_short("report_mols", std::ofstream::out | std::ofstream::app);
        std::ofstream out(output_filename);
        std::ofstream out_csv("report.csv", std::ofstream::out | std::ofstream::app);

        out_short << mol->get_file_name() << ":  ";
        for (int i = 0; i < nrpsMatchs.size(); ++i) {
            nrpsMatchs[i].print(out);
            nrpsMatchs[i].print_csv(out_csv);
            if (i < 3) {
                nrpsMatchs[i].print_short_prediction(out_short);
            }
        }
        out_short << "\n";

        out_short.close();
        out_csv.close();
        out.close();
    };
    delete(score);
}

std::string gen_filename(std::string ifile, std::string prefix) {
    size_t pos_sl = ifile.find_last_of('/');
    return (prefix + ifile.substr(pos_sl + 1) + ".match");
}

int main(int argc, char* argv[]) {
    logging::create_console_logger("");

    cxxopts::Options options = parse_options(argc, argv);

    std::string AA_file_name = options["amino_acids"].as<std::string>();
    std::string cfg_filename = options["config"].as<std::string>();
    Args args(cfg_filename);

    int start_from = std::max(0, options["start_from"].as<int>());
    aminoacid::AminoacidInfo::init(AA_file_name, "NRPSPREDICTOR2", args.aminoacid_info_default_logp);
    aminoacid::ModificationInfo::init(args.modification_cfg);
    aminoacid::MonomerInfo::init(args.monomer_cfg, args.monomer_logP_cfg, args.monomer_info_default_logp);

    std::cout << args.modification_cfg << "\n";
    for (int i = 0; i != aminoacid::ModificationInfo::MODIFICATION_CNT; ++i) {
        std::cout << aminoacid::ModificationInfo::NAMES[i] << " ";
        for (double x : aminoacid::ModificationInfo::COEFFICIENT[i]) {
            std::cout << x << " ";
        }
        std::cout << std::endl;
    }

    INFO("NRPs Matcher START");
    INFO("Saving predictions");
    std::vector<nrpsprediction::BgcPrediction> preds = save_predictions(options["predictions"].as<std::string>());
    INFO("Saving NRPs structures");
    std::vector<std::shared_ptr<nrp::NRP>> mols = load_nrps_from_monomeric_info(options["structures"].as<std::string>());

    if (start_from == 0) {
        std::ofstream out_csv("report.csv");
        out_csv << "score,peptide,nrp len,match cnt,all matched,mol id,prediction id\n";
        out_csv.close();
    }

    INFO("Processing matching for NRPs structurs")
    INFO("Start from: " << start_from)
    unsigned nthreads = args.threads;

    nerpa_set_omp_threads(nthreads);

    INFO("THREADS #" << nthreads);
#   pragma omp parallel for
    for (int i = start_from; i < mols.size(); ++i) {
        INFO("NRP structure #" << i)
        std::string output_filename = gen_filename(mols[i]->get_file_name(), "details_mols/");

        run_mol_predictions(preds, mols[i], output_filename, args);
    }

    return 0;
}