#include <iostream>
#include <fstream>
#include <NRP/NRPBuilder.h>
#include <algorithm>
#include <sstream>
#include "NRP/NRP.h"
#include "NRPsPrediction/NRPsPrediction.h"

const std::string MODE_PREDICTION_MOLS = "prediction_mols";
const std::string MODE_MOL_PREDICTIONS = "mol_predictions";

void run_prediction_mols_mode(char* argv[]) {
    std::string prediction_file(argv[2]);
    std::string path_to_nrps_file(argv[3]);

    std::ofstream out("nrpsMatch");
    std::ofstream out_short("report_predictions", std::ofstream::out | std::ofstream::app);

    out_short << prediction_file << ":  ";

    nrpsprediction::NRPsPrediction nrPsPrediction;
    nrPsPrediction.read_file(prediction_file);

    std::ifstream in_nrps_files(path_to_nrps_file);

    std::string cur_nrp_file;
    std::string cur_line;
    std::vector<nrp::NRP::Match> nrpsMatchs;
    std::vector<nrp::NRP*> nrpptr;
    while(getline(in_nrps_files, cur_line)) {
        std::stringstream ss(cur_line);
        ss >> cur_nrp_file;
        std::string extra_info;
        getline(ss, extra_info);
        nrp::NRP* nrp_from_fragment_graph = nrp::NRPBuilder::build(cur_nrp_file, extra_info);
        if (nrp_from_fragment_graph == nullptr) {
            continue;
        }

        nrpptr.push_back(nrp_from_fragment_graph);
        nrpsMatchs.push_back(nrp_from_fragment_graph->isCover(nrPsPrediction));
    }

    std::sort(nrpsMatchs.begin(), nrpsMatchs.end());
    for (int i = 0; i < nrpsMatchs.size(); ++i) {
        nrpsMatchs[i].print(out);
        if (i < 3) {
            nrpsMatchs[i].print_short(out_short);
        }
    }
    out_short << "\n";
    for (int i = 0; i < nrpptr.size(); ++i) {
        delete nrpptr[i];
    }

    out_short.close();
    out.close();
}

void run_mol_predictions_mode(char* argv[]) {
    std::string path_to_predictions_file(argv[2]);
    std::string nrp_file(argv[3]);

    std::ofstream out("nrpsMatch");
    std::ofstream out_short("report_mols", std::ofstream::out | std::ofstream::app);

    nrp::NRP* nrp_from_fragment_graph = nrp::NRPBuilder::build(nrp_file, "");
    if (nrp_from_fragment_graph == nullptr) {
        return;
    }

    out_short << nrp_file << ":  ";

    std::ifstream in_predictions_files(path_to_predictions_file);

    std::string cur_prediction_file;
    std::string cur_line;
    std::vector<nrp::NRP::Match> nrpsMatchs;
    while(getline(in_predictions_files, cur_line)) {
        std::stringstream ss(cur_line);
        ss >> cur_prediction_file;

        nrpsprediction::NRPsPrediction nrPsPrediction;
        nrPsPrediction.read_file(cur_prediction_file);

        nrpsMatchs.push_back(nrp_from_fragment_graph->isCover(nrPsPrediction));
    }

    std::sort(nrpsMatchs.begin(), nrpsMatchs.end());
    for (int i = 0; i < nrpsMatchs.size(); ++i) {
        nrpsMatchs[i].print(out);
        if (i < 3) {
            nrpsMatchs[i].print_short_prediction(out_short);
        }
    }

    out_short << "\n";
    delete nrp_from_fragment_graph;

    out_short.close();
    out.close();
}

int main(int argc, char* argv[]) {
    if (argv[1] == MODE_PREDICTION_MOLS) {
        run_prediction_mols_mode(argv);
    } else if (argv[1] == MODE_MOL_PREDICTIONS) {
        run_mol_predictions_mode(argv);
    }
    return 0;
}