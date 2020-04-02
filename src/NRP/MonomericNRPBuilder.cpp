//
// Created by tag on 31/03/2020.
//

#include "Aminoacid/MonomerInfo.h"
#include "MonomericNRPBuilder.h"
#include "NRPCycle.h"
#include "NRPLine.h"
#include "NRPtail.h"
#include "NRP.h"

std::shared_ptr<nrp::NRP> nrp::MonomericNRPBuilder::build(std::string nrp_id, std::string extra) {
    std::vector<aminoacid::Aminoacid> aminoacids;
    std::stringstream ss(extra);
    std::string graph, extra_info;
    ss >> graph;
    getline(ss, extra_info);

    size_t pos = graph.find(';');
    if (pos == std::string::npos) return nullptr;

    std::string monomers = graph.substr(0, pos);
    std::string bonds = graph.substr(pos + 1);

    ss = std::stringstream(monomers);
    std::string code;
    while (std::getline(ss, code, ',')) {
        aminoacids.emplace_back(aminoacid::MonomerInfo::getAAByCode(code));
    }

    std::vector<std::vector<int> > g(aminoacids.size());
    std::vector<std::vector<int> > gr(aminoacids.size());
    std::vector<std::string> strformula(aminoacids.size());

    ss = std::stringstream(bonds);
    std::string bond;
    std::string tmp;
    int b, e;
    while (std::getline(ss, bond, ';')) {
        std::stringstream ss_(bond);
        std::getline(ss_, tmp, ',');
        std::istringstream(tmp) >> b;
        std::getline(ss_, tmp, ',');
        std::istringstream(tmp) >> e;

        g[b].push_back(e);
        gr[e].push_back(b);
    }

    if (!isConnected(g, gr)) {
        return nullptr;
    }

    if (isCycle(g, gr)) {
        std::vector<int> pos = parseCycle(g, gr);
        std::vector<aminoacid::Aminoacid> resaacid = aminoacids_by_pos(aminoacids, pos);

        return std::make_shared<nrp::NRPCycle>(nrp_id, strformula, resaacid, pos, graph, extra_info);
    } else if (isLine(g, gr)) {
        std::vector<int> pos = parseLine(g, gr);
        std::vector<aminoacid::Aminoacid> resaacid = aminoacids_by_pos(aminoacids, pos);

        return std::make_shared<nrp::NRPLine>(nrp_id, strformula, resaacid, pos, graph, extra_info);
    } else if (isTail(g, gr)) {
        std::vector<int> pos_tail, pos_cycle;
        parseTail(g, gr, pos_tail, pos_cycle);

        std::vector<int> pos1 = pos_tail, pos2 = pos_tail;
        for (int i = 0; i < pos_cycle.size(); ++i) {
            pos1.push_back(pos_cycle[i]);
        }

        for (int i = pos_cycle.size() - 1; i >= 0; --i) {
            pos2.push_back(pos_cycle[i]);
        }

        std::vector<aminoacid::Aminoacid> resaacid1 = aminoacids_by_pos(aminoacids, pos1);
        std::vector<aminoacid::Aminoacid> resaacid2 = aminoacids_by_pos(aminoacids, pos2);

        std::shared_ptr<NRP> ver1 = std::make_shared<NRPLine>(nrp_id, strformula, resaacid1, pos1, graph, extra_info),
                ver2 = std::make_shared<NRPLine>(nrp_id, strformula, resaacid2, pos2, graph, extra_info);

        return std::make_shared<nrp::NRPtail>(ver1, ver2);
    }

    return nullptr;
}
