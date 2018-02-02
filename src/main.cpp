#include <iostream>
#include <fstream>
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
        nrp::NRP nrp_from_fragment_graph;
        nrp_from_fragment_graph.parse_fragment_graph(cur_nrp_file);

        nrpsmatch::NRPsMatch nrPsMatch = nrPsPrediction.isCover(nrp_from_fragment_graph);
        nrPsMatch.print(out);
    }

    out.close();
    return 0;
}