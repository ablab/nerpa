#include <iostream>
#include <fstream>
#include <NRP/NRPBuilder.h>
#include <algorithm>
#include <sstream>
#include <NRPGenerator/NRPGenerator.h>
#include <NormalizedMatch/NormalizedMatch.h>
#include <NRPGenerator/NRPGeneratorTriplet.h>
#include "NRP/NRP.h"
#include "NRPsPrediction/NRPsPrediction.h"

const int MIN_SCROE = 2;

std::vector<nrpsprediction::NRPsPrediction>  save_predictions(char* file_name) {
    std::vector<nrpsprediction::NRPsPrediction> preds;
    std::ifstream in_predictions_files(file_name);

    std::string cur_prediction_file;
    std::string cur_line;

    while(getline(in_predictions_files, cur_line)) {
        //std::stringstream ss(cur_line);
        //ss >> cur_prediction_file;

        nrpsprediction::NRPsPrediction nrPsPrediction;
        nrPsPrediction.read_file(cur_line);

        preds.push_back(nrPsPrediction);
    }

    return preds;
}

std::vector<nrp::NRP*> save_mols(char* file_name) {
    std::vector<nrp::NRP*> mols;

    std::ifstream in_nrps_files(file_name);
    std::string cur_nrp_file;
    std::string cur_line;

    while(getline(in_nrps_files, cur_line)) {
        std::stringstream ss(cur_line);
        ss >> cur_nrp_file;
        std::string extra_info;
        getline(ss, extra_info);
        nrp::NRP* nrp_from_fragment_graph = nrp::NRPBuilder::build(cur_nrp_file, extra_info);
        if (nrp_from_fragment_graph == nullptr) {
            continue;
        }
        mols.push_back(nrp_from_fragment_graph);
    }

    return mols;
}

void run_prediction_mols(nrpsprediction::NRPsPrediction pred, std::vector<nrp::NRP*> mols,
                         nrp_generator::NRPGenerator* nrpGenerator, std::string output_filename) {
    if (pred.getNrpsParts().size() == 0) return;
    std::ofstream out(output_filename);
    std::ofstream out_short("report_predictions", std::ofstream::out | std::ofstream::app);

    std::vector<normalized_match::NormalizedMatch> nrpsMatchs;
    for (int i = 0; i < mols.size(); ++i) {
        nrp::NRP::Match match = mols[i]->isCover(pred);
        if (match.score() >= MIN_SCROE) {
            nrpsMatchs.push_back(normalized_match::NormalizedMatch(match, nrpGenerator, pred, mols[i]));
        }
    }

    if (nrpsMatchs.size() > 0) {
        out_short << pred.getNrpsParts()[0].get_file_name() << ":  ";
    }
    std::sort(nrpsMatchs.begin(), nrpsMatchs.end());
    for (int i = 0; i < nrpsMatchs.size(); ++i) {
        nrpsMatchs[i].print(out);
        if (i < 3) {
            nrpsMatchs[i].print_short(out_short);
        }
    }
    if (nrpsMatchs.size() > 0) {
        out_short << "\n";
    }

    out_short.close();
    out.close();
}

void run_mol_predictions(std::vector<nrpsprediction::NRPsPrediction> preds, nrp::NRP* mol,
                         nrp_generator::NRPGenerator* nrpGenerator, std::string output_filename) {
    std::ofstream out_short("report_mols", std::ofstream::out | std::ofstream::app);

    std::vector<normalized_match::NormalizedMatch> nrpsMatchs;
    for (int i = 0; i < preds.size(); ++i) {
        nrp::NRP::Match match = mol->isCover(preds[i]);
        if (match.score() >= MIN_SCROE) {
            nrpsMatchs.push_back(normalized_match::NormalizedMatch(match, nrpGenerator, preds[i], mol));
        }
    }

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
    std::vector<nrpsprediction::NRPsPrediction> preds = save_predictions(argv[1]);
    std::vector<nrp::NRP*> mols = save_mols(argv[2]);
    nrp_generator::NRPGenerator* nrpGenerator = new nrp_generator::NRPGeneratorTriplet(mols);
    for (int i = 0; i < preds.size(); ++i) {
        std::cerr << "pred " << i << "\n";
        if (preds[i].getNrpsParts().size() == 0) continue;
        std::string output_filename = gen_filename(preds[i].getNrpsParts()[0].get_file_name(), "details_predictions/");
        run_prediction_mols(preds[i], mols, nrpGenerator, output_filename);
    }
    std::ofstream out_csv("report.csv");
    out_csv << "score,normalize score,peptide,nrp len,match cnt,all matched,mol id,prediction id, p-value\n";
    out_csv.close();

    for (int i = 0; i < mols.size(); ++i) {
        std::cerr << "mol " << i << "\n";
        std::string output_filename = gen_filename(mols[i]->get_file_name(), "details_mols/");
        run_mol_predictions(preds, mols[i], nrpGenerator, output_filename);
    }

    delete nrpGenerator;
    return 0;
}