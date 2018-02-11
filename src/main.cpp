#include <iostream>
#include <fstream>
#include <NRP/NRPBuilder.h>
#include "NRP/NRP.h"
#include "NRPsPrediction/NRPsPrediction.h"

int main(int argc, char* argv[]) {
    std::string prediction_file(argv[1]);
    std::string path_to_nrps_file(argv[2]);

    std::ofstream out("nrpsMatch");

    nrpsprediction::NRPsPrediction nrPsPrediction;
    nrPsPrediction.read_file(prediction_file);

    std::ifstream in_nrps_files(path_to_nrps_file);

    std::string cur_nrp_file;
    while(getline(in_nrps_files, cur_nrp_file)) {
        std::cerr << cur_nrp_file << "\n";
        nrp::NRP* nrp_from_fragment_graph = nrp::NRPBuilder::build(cur_nrp_file);
        if (nrp_from_fragment_graph == nullptr) {
            continue;
        }
        nrp::NRP::Match nrPsMatch = nrp_from_fragment_graph->isCover(nrPsPrediction);
        nrPsMatch.print(out);
        delete nrp_from_fragment_graph;
    }

    std::cerr << "ok8\n";
    out.close();
    return 0;
}