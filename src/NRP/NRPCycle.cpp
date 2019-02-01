#include <iostream>
#include <Logger/logger.hpp>
#include "NRPCycle.h"

nrp::NRPCycle::NRPCycle(const std::string &file_name, const std::vector<std::string> &strformula,
                        const std::vector<aminoacid::Aminoacid> &aminoacids,
                        const std::vector<int> &position, const std::string &graph, const std::string& extra_info) : NRP(file_name, strformula,
                                                                                          aminoacids, position,
                                                                                          graph, extra_info) {
}

nrp::NRP::NRPType nrp::NRPCycle::getType() const {
    return NRP::cycle;
}

//TODO
std::vector<nrp::NRPLine> nrp::NRPCycle::getLines() const {
    ERROR("Get lines for NRP cycle. NOT implemented");
    return std::vector<nrp::NRPLine>();
}

bool nrp::NRPCycle::is_valid_seg(int l, int r, int stp) const {
    return true;
}
