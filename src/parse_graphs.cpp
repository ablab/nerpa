#include <iostream>
#include <fstream>
#include <NRP/NRPBuilder.h>
#include <algorithm>
#include <Logger/log_writers.hpp>

std::string get_file_name(std::string cur_line) {
    std::string res = "";
    for (int i = 0; i < cur_line.size() && cur_line[i] != '\t'; ++i) {
        res += cur_line[i];
    }
    return res;
}

std::vector<std::shared_ptr<nrp::NRP>> save_mols(char* file_name) {
    std::vector<std::shared_ptr<nrp::NRP>> mols;

    std::ifstream in_nrps_files(file_name);
    std::string cur_nrp_file;
    std::string cur_line;

    while(getline(in_nrps_files, cur_line)) {
        std::stringstream ss(cur_line);
        ss >> cur_nrp_file;
        std::string extra_info;
        getline(ss, extra_info);
        std::shared_ptr<nrp::NRP> nrp_from_fragment_graph = nrp::NRPBuilder::build(cur_nrp_file, extra_info);
        if (nrp_from_fragment_graph == nullptr) {
            std::cerr << cur_nrp_file << " " << -1 << " " << 0 << "\n";
            continue;
        }
        std::string typeName[3] = {"CYCLE", "LINE", "BRANCH"};

        std::cerr << cur_nrp_file << " " << typeName[nrp_from_fragment_graph->getType()] << " " << nrp_from_fragment_graph->getLen() << "\n";
        nrp_from_fragment_graph->print();

        mols.push_back(nrp_from_fragment_graph);
    }

    return mols;
}


int main(int argc, char* argv[]) {
    logging::create_console_logger("");

    aminoacid::AminoacidInfo::init(argv[2], "MINOWA");
    std::vector<std::shared_ptr<nrp::NRP>> mols = save_mols(argv[1]);

    return 0;
}