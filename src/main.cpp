#include <iostream>
#include <fstream>
#include <NRP/NRPBuilder.h>
#include <algorithm>
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
    std::vector<nrp::NRP::Match> nrpsMatchs;
    std::vector<nrp::NRP*> nrpptr;
    while(getline(in_nrps_files, cur_nrp_file)) {
        std::cerr << cur_nrp_file << "\n";
        nrp::NRP* nrp_from_fragment_graph = nrp::NRPBuilder::build(cur_nrp_file);
        if (nrp_from_fragment_graph == nullptr) {
            continue;
        }

        nrpptr.push_back(nrp_from_fragment_graph);
        nrpsMatchs.push_back(nrp_from_fragment_graph->isCover(nrPsPrediction));
    }

    std::sort(nrpsMatchs.begin(), nrpsMatchs.end());
    for (int i = 0; i < nrpsMatchs.size(); ++i) {
        nrpsMatchs[i].print(out);
    }

    for (int i = 0; i < nrpptr.size(); ++i) {
        delete nrpptr[i];
    }

    out.close();
    return 0;
}