#include <iostream>
#include <fstream>
#include <NRP/NRPBuilder.h>
#include <NRP/MonomericNRPBuilder.h>
#include <algorithm>
#include <sstream>
#include <cstring>
#include "NRP/NRP.h"
#include "NRPsPrediction/BgcPrediction.h"
#include <Logger/log_writers.hpp>
#include <NRPsPrediction/Builders/Nrpspredictor2Builder.h>
#include <Matcher/Score/Base/ScoreWithModification.h>
#include <NRPsPrediction/Builders/MinowaPredictionBuilder.h>
#include <Matcher/Score/Minowa/ScoreMinowa.h>
#include <NRPsPrediction/Builders/PrismPredictionBuilder.h>
#include <Matcher/Score/Prism/ScorePrism.h>
#include <Matcher/Score/Sandpuma/ScoreSandpuma.h>
#include <NRPsPrediction/Builders/SandpumaPredictionBuilder.h>
#include <Matcher/Score/Base/ScorePositionOnly.h>
#include <Matcher/Score/Minowa/ScoreMinowaScoreOnly.h>
#include <Matcher/Score/NrpsPredictor2/ScoreNRPsPredictor2Normalize.h>
#include <Matcher/Score/Minowa/ScoreMinowaPositionalCoefficient.h>
#include <Matcher/Score/Base/ScoreSingleUnit.h>
#include <Matcher/Score/Base/ScoreOpenContinueGap.h>
#include <Matcher/Score/Base/ScoreNormalize.h>
#include <Matcher/Score/Base/ScoreForUntrustedPred.h>
#include <ArgParse/Args.h>
#include <Matcher/SingleUnitMatcher.h>
#include <Aminoacid/ModificationInfo.h>
#include <omp.h>
#include <Matcher/OrderedGenesMatcher.h>
#include <Matcher/Score/OrderedGenes/OrderedGenesScoreBase.h>
#include <Aminoacid/MonomerInfo.h>
#include "Matcher/Matcher.h"
#include "Matcher/InDelMatcher.h"

const double MIN_SCROE = 0.05;
const double MIN_EXPLAIN_PART = 0;//0.15;

void getPredictor(std::string predictor_name, nrpsprediction::PredictionBuilderBase*& predictionBuilder) {
    if (predictor_name == "MINOWA") {
        predictionBuilder = new nrpsprediction::MinowaPredictionBuilder();
    } else if (predictor_name == "PRISM") {
        predictionBuilder = new nrpsprediction::PrismPredictionBuilder();
    } else if (predictor_name == "SANDPUMA") {
        predictionBuilder = new nrpsprediction::SandpumaPredictionBuilder();
    } else if (predictor_name == "NRPSPREDICTOR2"){
        predictionBuilder = new nrpsprediction::Nrpspredictor2Builder();
    } else {
        ERROR("Unknown predictor " + predictor_name);
    }
}

std::string get_file_name(std::string cur_line) {
    std::string res = "";
    for (int i = 0; i < cur_line.size() && cur_line[i] != '\t'; ++i) {
        res += cur_line[i];
    }
    return res;
}

std::vector<nrpsprediction::BgcPrediction>  save_predictions(char* file_name, std::string predictor_name) {
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
        getPredictor(predictor_name, nrPsPredictionBuilder);
        nrPsPredictionBuilder->read_file(info_file_name);

        preds.push_back(nrPsPredictionBuilder->getPrediction());
        INFO("Parts in prediction: " << preds.back().getOrfs().size());
        delete(nrPsPredictionBuilder);
    }

    return preds;
}

std::vector<std::shared_ptr<nrp::NRP>> save_mols(char* file_name) {
    std::vector<std::shared_ptr<nrp::NRP>> mols;

    std::ifstream in_nrps_files(file_name);
    std::ofstream out_csv("structure_details.csv", std::ofstream::out);
    out_csv << "Accession\tSTRUCTURE\tVERTEX\tFORMULA\tAA\n";

    std::string cur_nrp_file;
    std::string cur_line;

    while(getline(in_nrps_files, cur_line)) {
        INFO(cur_line);
        std::stringstream ss(cur_line);
        ss >> cur_nrp_file;
        std::string extra_info;
        getline(ss, extra_info);
        std::shared_ptr<nrp::NRP> nrp_from_fragment_graph = nrp::NRPBuilder::build(cur_nrp_file, extra_info);
        if (nrp_from_fragment_graph == nullptr) {
            continue;
        }
        nrp_from_fragment_graph->print();
        for (int i = 0; i < nrp_from_fragment_graph->getFullLen(); ++i) {
            out_csv << cur_line << "\t" << nrp_from_fragment_graph->structure_to_string() << "\t" << i
            << "\t" << nrp_from_fragment_graph->getFormula(i) << "\t"
            << nrp_from_fragment_graph->getAminoacid(i).get_possible_name() << "\n";
        }
        mols.push_back(nrp_from_fragment_graph);
    }

    out_csv.close();
    return mols;
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
    if (args.predictor_name == "MINOWA") {
        score = new ScoreMinowa(args.mismatch);
    } else if (args.predictor_name == "PRISM") {
        score = new ScorePrism(args.mismatch);
    } else if (args.predictor_name == "SANDPUMA") {
        score = new ScoreSandpuma(args.mismatch);
    } else {
        score = new Score(args.mismatch);
    }
    score = new OrderedGenesScoreBase(std::unique_ptr<Score>(std::move(score)), args.skip_segment, args.insertion, args.deletion, args.mismatch, args.open_gap, args.continue_gap);

    if (args.modification) {
        score = new ScoreWithModification(std::unique_ptr<Score>(std::move(score)));
    }
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
    int pos_sl = -1;
    int pos_pt = -1;
    for (int i = 0; i < ifile.size(); ++i) {
        if (ifile[i] == '/') {
            pos_sl = i;
        }
        if (ifile[i] == '.') {
            pos_pt = i;
        }
    }

    return (prefix + ifile.substr(pos_sl + 1,  pos_pt - pos_sl - 1));
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
    aminoacid::AminoacidInfo::init(AA_file_name, args.predictor_name);
    aminoacid::ModificationInfo::init(args.modification_cfg);
    aminoacid::ModificationInfo::init_AAMod(args.AAmod_cfg);
    aminoacid::MonomerInfo::init(args.monomer_cfg, args.monomer_logP_cfg);

    INFO("NRPs Matcher START");
    INFO("Saving predictions");
    std::vector<nrpsprediction::BgcPrediction> preds = save_predictions(argv[1], args.predictor_name);
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
#   pragma omp parallel for
    for (int i = start_from; i < mols.size(); ++i) {
        INFO("NRP structure #" << i)
        std::string output_filename = gen_filename(mols[i]->get_file_name(), "details_mols/");

        run_mol_predictions(preds, mols[i], output_filename, args);
    }

    return 0;
}