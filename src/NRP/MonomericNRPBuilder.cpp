//
// Created by tag on 31/03/2020.
//

#include <Logger/logger.hpp>
#include "Aminoacid/MonomerInfo.h"
#include "MonomericNRPBuilder.h"
#include "NRPCycle.h"
#include "NRPLine.h"
#include "NRPBranch.h"
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

    std::vector<std::string> strnodes;
    ss = std::stringstream(monomers);
    std::string code;
    size_t node_idx = 0;
    while (std::getline(ss, code, ',')) {
        aminoacid::Aminoacid aa;
        if (code[0] == '*') {
            std::string code_ = code.substr(1);
            aa = aminoacid::MonomerInfo::getAAByCode(code_);
            // TODO: consider introducing a special modification type
            aminoacid::Modification mod(aminoacid::ModificationInfo::MODIFICATION_CNT - 1);
            aa.addModification(mod);
        } else {
            aa = aminoacid::MonomerInfo::getAAByCode(code);
        }

        aminoacid::Aminoacid::Configuation cur_config = aminoacid::Aminoacid::NA;
        for (int i = 0 ; i < code.size(); ++i) {
            if (code[i] == '@' && code[i + 1] == 'D') {
                cur_config = aminoacid::Aminoacid::D;
            }

            if (code[i] == '@' && code[i + 1] == 'L') {
                cur_config = aminoacid::Aminoacid::L;
            }

            if (code[i] == '@' && code[i + 1] == 'P') {
                aminoacid::Modification mod(aminoacid::ModificationInfo::getIdByNameId("pkhybrid"));
                aa.addModification(mod);
            }
        }

        aa.setConfiguration(cur_config);
        aminoacids.push_back(aa);

        std::stringstream ssnode;
        ssnode << node_idx++ << " " << code << " " << 0;
        strnodes.push_back(ssnode.str());
    }

    std::vector<std::vector<int> > g(aminoacids.size());
    std::vector<std::vector<int> > gr(aminoacids.size());

    ss = std::stringstream(bonds);
    std::string bond;
    std::string tmp;
    int b, e;
    std::vector<std::pair<int, int>> ocon;
    while (std::getline(ss, bond, ';')) {
        std::stringstream ss_(bond);

        if (bond.empty()) break;
        
        std::getline(ss_, tmp, ',');
        std::istringstream(tmp) >> b;
        std::getline(ss_, tmp, ',');
        bool is_ocon = false;
        if (tmp == "o") {
            is_ocon = true;
            std::getline(ss_, tmp, ',');
        }
        std::istringstream(tmp) >> e;

        g[b].push_back(e);
        gr[e].push_back(b);
        if (is_ocon) {
            ocon.emplace_back(std::make_pair(b, e));
        }
    }

    std::string strgraph = graphToString(g);

    if (!isConnected(g, gr)) {
        return nullptr;
    }

    if (isCycle(g, gr)) {
        std::vector<int> pos = parseCycle(g, gr, ocon);
        std::vector<aminoacid::Aminoacid> resaacid = aminoacids_by_pos(aminoacids, pos);
        if (ocon.empty()) {
            return std::make_shared<nrp::NRPCycle>(nrp_id, strnodes, resaacid, pos, strgraph, extra_info);
        } else {
            return std::make_shared<nrp::NRPLine>(nrp_id, strnodes, resaacid, pos, strgraph, extra_info);
        }
    } else if (isLine(g, gr)) {
        std::vector<int> pos = parseLine(g, gr);
        std::vector<aminoacid::Aminoacid> resaacid = aminoacids_by_pos(aminoacids, pos);

        return std::make_shared<nrp::NRPLine>(nrp_id, strnodes, resaacid, pos, strgraph, extra_info);
    } else if (isTail(g, gr)) {
        std::vector<int> pos_tail, pos_cycle;
        parseTail(g, gr, pos_tail, pos_cycle, ocon);

        std::vector<int> pos1 = pos_tail, pos2 = pos_tail;
        for (int i = 0; i < pos_cycle.size(); ++i) {
            pos1.push_back(pos_cycle[i]);
        }

        bool is_one_pos = false;
        for (int i = 0; i < ocon.size(); ++i) {
            if (ocon[i].first == pos_tail.back() && ocon[i].second == pos_cycle.back()) {
                is_one_pos = true;
            }

            if (ocon[i].first == pos_cycle.back() && ocon[i].second == pos_tail.back()) {
                is_one_pos = true;
            }
        }
        std::vector<aminoacid::Aminoacid> resaacid1 = aminoacids_by_pos(aminoacids, pos1);

        if (is_one_pos) {
            return std::make_shared<NRPLine>(nrp_id, strnodes, resaacid1, pos1, strgraph, extra_info);
        }

        for (int i = pos_cycle.size() - 1; i >= 0; --i) {
            pos2.push_back(pos_cycle[i]);
        }

        std::vector<aminoacid::Aminoacid> resaacid2 = aminoacids_by_pos(aminoacids, pos2);

        std::shared_ptr<NRP> ver1 = std::make_shared<NRPLine>(nrp_id, strnodes, resaacid1, pos1, strgraph, extra_info),
                ver2 = std::make_shared<NRPLine>(nrp_id, strnodes, resaacid2, pos2, strgraph, extra_info);

        return std::make_shared<nrp::NRPBranch>(ver1, ver2, pos_tail.size());
    }

    return nullptr;
}
