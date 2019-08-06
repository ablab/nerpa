#include <iostream>
#include <fstream>
#include <NRP/NRPBuilder.h>
#include <algorithm>
#include <sstream>
#include <cstring>
#include "NRP/NRP.h"
#include "NRPsPrediction/NRPsPrediction.h"
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
#include <ArgParse/Args.h>
#include <Matcher/SingleUnitMatcher.h>
#include "Matcher/Matcher.h"
#include "Matcher/InDelMatcher.h"

const double MIN_SCROE = 0.002;
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

std::vector<nrpsprediction::NRPsPrediction>  save_predictions(char* file_name, std::string predictor_name) {
    std::vector<nrpsprediction::NRPsPrediction> preds;
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
        INFO("Parts in prediction: " << preds.back().getNrpsParts().size());
        delete(nrPsPredictionBuilder);
    }

    return preds;
}

std::vector<std::shared_ptr<nrp::NRP>> save_mols(char* file_name) {
    std::vector<std::shared_ptr<nrp::NRP>> mols;

    std::ifstream in_nrps_files(file_name);
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
        mols.push_back(nrp_from_fragment_graph);
    }

    return mols;
}

void getScoreFunction(Args args, matcher::Score*& score) {
    using namespace matcher;
    if (args.predictor_name == "MINOWA") {
        score = new ScoreMinowa;
    } else if (args.predictor_name == "PRISM") {
        score = new ScorePrism;
    } else if (args.predictor_name == "SANDPUMA") {
        score = new ScoreSandpuma;
    } else {
        score = new Score;
    }

    score = new ScoreNormalize(std::unique_ptr<Score>(std::move(score)));
    score = new ScoreOpenContinueGap(args.open_gap, args.continue_gap,
                                     std::unique_ptr<Score>(std::move(score)));
    if (args.single_match) {
        score = new ScoreSingleUnit(args.single_match_coeff, std::unique_ptr<Score>(std::move(score)));
    }
    if (args.modification) {
        score = new ScoreWithModification(std::unique_ptr<Score>(std::move(score)));
    }
}

matcher::MatcherBase* getMatcher(Args args) {
    matcher::MatcherBase* matcherBase;
    if (args.single_match) {
        matcherBase = new matcher::SingleUnitMatcher();
    } else {
        matcherBase = new matcher::Matcher();
    }
    if (args.deletion || args.insertion) {
        return new matcher::InDelMatcher(matcherBase, args.insertion, args.deletion);
    } else {
        return matcherBase;
    }
}


void run_mol_predictions(std::vector<nrpsprediction::NRPsPrediction> preds, std::shared_ptr<nrp::NRP> mol, std::string output_filename,
                         Args args) {
    std::ofstream out_short("report_mols", std::ofstream::out | std::ofstream::app);
    matcher::Score* score;
    getScoreFunction(args, score);
    std::vector<matcher::MatcherBase::Match> nrpsMatchs;
    for (int i = 0; i < preds.size(); ++i) {
        if (preds[i].getNrpsParts().size() == 0) continue;
        //std::cerr << mol->get_file_name() << " " << preds[i].getNrpsParts()[0].get_file_name() << "\n";
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
    std::ofstream out(output_filename);
    out_short << mol->get_file_name() << ":  ";

    std::ofstream out_csv("report.csv", std::ofstream::out | std::ofstream::app);
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

    INFO("NRPs Matcher START");
    INFO("Saving predictions");
    std::vector<nrpsprediction::NRPsPrediction> preds = save_predictions(argv[1], args.predictor_name);
    INFO("Saving NRPs structures");
    std::vector<std::shared_ptr<nrp::NRP>> mols = save_mols(argv[2]);

    if (start_from == 0) {
        std::ofstream out_csv("report.csv");
        out_csv << "score,peptide,nrp len,match cnt,all matched,mol id,prediction id\n";
        out_csv.close();
    }

    INFO("Processing matching for NRPs structurs")
    INFO("Start from: " << start_from)
    for (int i = start_from; i < mols.size(); ++i) {
        INFO("NRP structure #" << i)
        std::string output_filename = gen_filename(mols[i]->get_file_name(), "details_mols/");

        run_mol_predictions(preds, mols[i], output_filename, args);
    }

    return 0;
}