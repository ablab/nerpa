#include <iostream>
#include <fstream>
#include <NRP/NRPBuilder.h>
#include "NRP/NRP.h"
#include "NRPsPrediction/NRPsPrediction.h"

int main(int argc, char* argv[]) {
    std::cerr << "ok0\n";
    std::string prediction_file(argv[1]);
    std::string path_to_nrps_file(argv[2]);

    std::cerr << "ok1\n";
    std::ofstream out("nrpsMatch");

    std::cerr << "ok2\n";
    nrpsprediction::NRPsPrediction nrPsPrediction;
    nrPsPrediction.read_file(prediction_file);

    std::cerr << "ok3\n";
    std::ifstream in_nrps_files(path_to_nrps_file);

    std::cerr << "ok4\n";
    std::string cur_nrp_file;
    while(getline(in_nrps_files, cur_nrp_file)) {
        std::cerr << "ok5\n";
        nrp::NRP* nrp_from_fragment_graph = nrp::NRPBuilder::build(cur_nrp_file);
        std::cerr << "ok6\n";
        nrp::NRP::Match nrPsMatch = nrp_from_fragment_graph->isCover(nrPsPrediction);
        std::cerr << "ok7\n";
        nrPsMatch.print(out);
        delete nrp_from_fragment_graph;
    }

    std::cerr << "ok8\n";
    out.close();
    return 0;
}