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
#include <omp.h>
#include <Matcher/OrderedGenesMatcher.h>
#include <Aminoacid/MonomerInfo.h>

const double MIN_SCROE = 0.05;
const double MIN_EXPLAIN_PART = 0;//0.15;

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

std::vector<nrpsprediction::BgcPrediction>  save_predictions(char* file_name) {
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

std::vector<std::shared_ptr<nrp::NRP>> load_nrps_from_monomeric_info(char* file_name) {
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
        std::shared_ptr<nrp::NRP> nrp_from_fragment_graph = nrp::MonomericNRPBuilder::build(cur_id, extra);
        if (nrp_from_fragment_graph == nullptr) {
            continue;
        }
        nrp_from_fragment_graph->print();
        nrps.push_back(nrp_from_fragment_graph);
    }
    return nrps;
}

void getScoreFunction(Args args, matcher::Score*& score) {
    using namespace matcher;
    score = new Score(args.insertion, args.deletion);
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

        if (match.score() >= MIN_SCROE &&
                (double)match.getCntMatch()/preds[i].getSumPredictionLen() >= MIN_EXPLAIN_PART) {
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

    std::string AA_file_name = argv[3];
    std::string cfg_filename = argv[4];
    Args args(cfg_filename);

    int start_from = 0;
    if (argc > 5) {
        std::stringstream ss(argv[5]);
        ss >> start_from;
    }
    aminoacid::AminoacidInfo::init(AA_file_name, "NRPSPREDICTOR2");
    aminoacid::ModificationInfo::init(args.modification_cfg);
    aminoacid::MonomerInfo::init(args.monomer_cfg, args.monomer_logP_cfg);

    INFO("NRPs Matcher START");
    INFO("Saving predictions");
    std::vector<nrpsprediction::BgcPrediction> preds = save_predictions(argv[1]);
    INFO("Saving NRPs structures");
    std::vector<std::shared_ptr<nrp::NRP>> mols = load_nrps_from_monomeric_info(argv[2]);

    if (start_from == 0) {
        std::ofstream out_csv("report.csv");
        out_csv << "score,peptide,nrp len,match cnt,all matched,mol id,prediction id\n";
        out_csv.close();
    }

    INFO("Processing matching for NRPs structurs")
    INFO("Start from: " << start_from)
    unsigned nthreads = args.threads;

    omp_set_dynamic(0);
    omp_set_num_threads(nthreads);

    INFO("THREADS #" << nthreads);
//#   pragma omp parallel for
    for (int i = start_from; i < mols.size(); ++i) {
        INFO("NRP structure #" << i)
        std::string output_filename = gen_filename(mols[i]->get_file_name(), "details_mols/");

        run_mol_predictions(preds, mols[i], output_filename, args);
    }

    return 0;
}