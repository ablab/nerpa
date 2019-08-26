#include <Logger/logger.hpp>
#include "NRPLine.h"

nrp::NRPLine::NRPLine(const std::string &file_name, const std::vector<std::string> &strformula,
                      const std::vector<aminoacid::Aminoacid> &aminoacids, const std::vector<int> &position,
                      const std::string &graph, const std::string& extra_info) : NRP(file_name, strformula, aminoacids, position, graph, extra_info) {}

nrp::NRP::NRPType nrp::NRPLine::getType() const {
    return NRP::line;
}

//TODO
std::vector<std::shared_ptr<nrp::NRP>> nrp::NRPLine::getLines() const {
    ERROR("Get lines for NRP line. NOT implemented");
    return std::vector<std::shared_ptr<NRP>>();
}

bool nrp::NRPLine::is_valid_seg(int l, int r, int stp) const {
    return (l <= r && stp > 0) || (r <= l && stp < 0);
}

nrp::NRPLine::NRPLine(const nrp::NRP &refNrp) : NRP(refNrp) {
}
