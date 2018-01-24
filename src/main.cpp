#include <iostream>
#include "NRP/NRP.h"
#include "NRPsPrediction/NRPsPrediction.h"

int main() {
    nrp::NRP nrp_from_fragment_graph;
    nrp_from_fragment_graph.parse_fragment_graph("fragmented_graph_2");
    nrpsprediction::NRPsPrediction nrPsPrediction;
    nrPsPrediction.read_file("ctg1_nrpspredictor2_codes.txt");
    std::vector<nrpsprediction::NRPsPart> nrpparts = nrPsPrediction.getNrpsParts();
    nrp_from_fragment_graph.print();
    std::cout << nrpparts.size() << "\n";
    for (int i = 0; i < (int)nrpparts.size(); ++i) {
        std::cout << nrp_from_fragment_graph.containNRPsPart(nrpparts[i]) << "\n";
    }
    return 0;
}