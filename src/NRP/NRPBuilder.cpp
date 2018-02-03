#include <fstream>
#include <sstream>
#include <iostream>
#include <cassert>
#include "NRPBuilder.h"
#include "NRPCycle.h"
#include "NRPLine.h"
#include "NRPtail.h"

nrp::NRP* nrp::NRPBuilder::build(std::string fragment_graph) {
    std::ifstream in(fragment_graph);
    std::string s;
    std::vector<std::string> strformula;
    std::vector<aminoacid::Aminoacids::Aminoacid> aminoacids;

    int line_cnt;
    in >> s >> s >> s >> s >> line_cnt;

    for (int i = 0; i < line_cnt; ++i) {
        int id;
        std::string formula;
        double mass;
        in >> id >> formula >> mass;
        std::stringstream ss;
        ss << id << " " << formula << " " << mass;
        strformula.push_back(ss.str());
        aminoacids.push_back(aminoacid::Aminoacids::get_aminoacid_from_formula(formula));
    }

    std::vector<std::vector<int> > g(aminoacids.size());
    std::vector<std::vector<int> > gr(aminoacids.size());

    std::string graph;
    std::string tmp;
    while(getline(in, tmp)) {
        if (tmp.size() > 0) {
            graph += tmp;
            graph += "\n";

            if (tmp[0] != 'n') {
                std::stringstream ss(tmp);
                int b, e;
                ss >> b >> tmp >> e;
                g[b].push_back(e);
                gr[e].push_back(b);
            }
        }
    }

    in.close();

    if (isCycle(g, gr)) {
        std::vector<int> pos = parseCycle(g, gr);
        std::vector<aminoacid::Aminoacids::Aminoacid> resaacid = aminoacids_by_pos(aminoacids, pos);

        return new nrp::NRPCycle(fragment_graph, strformula, resaacid, pos, graph);
    } else if (isLine(g, gr)) {
        std::vector<int> pos = parseLine(g, gr);
        std::vector<aminoacid::Aminoacids::Aminoacid> resaacid = aminoacids_by_pos(aminoacids, pos);

        return new NRPLine(fragment_graph, strformula, resaacid, pos, graph);
    } else if (isTail(g, gr)) {
        std::vector<int> pos_tail, pos_cycle;
        parseTail(g, gr, pos_tail, pos_cycle);

        std::vector<int> pos1 = pos_tail, pos2 = pos_tail;
        for (int i = 0; i < pos_cycle.size(); ++i) {
            pos1.push_back(pos_cycle[i]);
        }

        for (int i = pos_cycle.size() - 1; i >= 0; ++i) {
            pos2.push_back(pos_cycle[i]);
        }

        std::vector<aminoacid::Aminoacids::Aminoacid> resaacid1 = aminoacids_by_pos(aminoacids, pos1);
        std::vector<aminoacid::Aminoacids::Aminoacid> resaacid2 = aminoacids_by_pos(aminoacids, pos2);

        NRPLine ver1(fragment_graph, strformula, resaacid1, pos1, graph), ver2(fragment_graph, strformula, resaacid2, pos2, graph);

        return new NRPtail(ver1, ver2);
    }

    assert(false);
}

std::vector<aminoacid::Aminoacids::Aminoacid>
nrp::NRPBuilder::aminoacids_by_pos(const std::vector<aminoacid::Aminoacids::Aminoacid> &aminoacids,
                              const std::vector<int> &pos) {
    std::vector<aminoacid::Aminoacids::Aminoacid> resaacid;
    for (int i = 0; i < pos.size(); ++i) {
            resaacid.push_back(aminoacids[pos[i]]);
        }
    return resaacid;
}

bool nrp::NRPBuilder::isCycle(std::vector<std::vector<int>> &g, std::vector<std::vector<int>> &gr) {
    for (int i = 0; i < (int)g.size(); ++i) {
        if (g[i].size() != 1) {
            return false;
        }
    }

    for (int i = 0; i < (int)g.size(); ++i) {
        if (gr[i].size() != 1) {
            return false;
        }
    }

    return true;
}

bool nrp::NRPBuilder::isLine(std::vector<std::vector<int>> &g, std::vector<std::vector<int>> &gr) {
    int cnt01 = 0, cnt02 = 0;
    for (int i = 0; i < (int)g.size(); ++i) {
        if (g[i].size() == 0) {
            cnt01++;
        }
    }

    for (int i = 0; i < (int)gr.size(); ++i) {
        if (gr[i].size() == 0) {
            cnt02++;
        }
    }

    return (cnt01 == 1 && cnt02 == 1);
}

bool nrp::NRPBuilder::isTail(std::vector<std::vector<int>> &g, std::vector<std::vector<int>> &gr) {
    int cnt01 = 0, cnt02 = 0;
    for (int i = 0; i < (int)g.size(); ++i) {
        if (g[i].size() == 0) {
            cnt01++;
        }
    }

    for (int i = 0; i < (int)gr.size(); ++i) {
        if (gr[i].size() == 0) {
            cnt02++;
        }
    }

    return (cnt01 + cnt02 == 1);
}

std::vector<int> nrp::NRPBuilder::parseCycle(std::vector<std::vector<int>> &g, std::vector<std::vector<int>> &gr) {
    std::vector<int> pos;
    pos.push_back(0);
    int cur = 0;
    while (pos.size() < g.size()) {
        cur = g[cur][0];
        pos.push_back(cur);
    }

    return pos;
}

std::vector<int> nrp::NRPBuilder::parseLine(std::vector<std::vector<int>> &g, std::vector<std::vector<int>> &gr) {
    int strp1 = -1;
    for (int i = 0; i < (int) g.size(); ++i) {
        if (g[i].size() == 0) {
            strp1 = i;
        }
    }
    std::vector<int> pos;
    pos.push_back(strp1);
    int cur = strp1;
    while (pos.size() < g.size()) {
        cur = gr[cur][0];
        pos.push_back(cur);
    }

    return pos;
}

void
nrp::NRPBuilder::parseTail(std::vector<std::vector<int>> &g, std::vector<std::vector<int>> &gr, std::vector<int> &tail,
                           std::vector<int> &cycle) {

    int strp1 = -1, strp2 = -1;
    for (int i = 0; i < (int) g.size(); ++i) {
        if (g[i].size() == 0) {
            strp1 = i;
        }
    }

    for (int i = 0; i < (int) gr.size(); ++i) {
        if (gr[i].size() == 0) {
            strp2 = i;
        }
    }

    if (strp1 != -1) {
        tail.push_back(strp1);
        int cur = strp1;
        while (gr[cur].size() == 1) {
            cur = gr[cur][0];
            tail.push_back(cur);
        }

        int nc = gr[cur][0];
        while (nc != cur) {
            cycle.push_back(nc);
            nc = gr[nc][0];
        }
    } else {
        tail.push_back(strp2);
        int cur = strp2;
        while (g[cur].size() == 1) {
            cur = g[cur][0];
            tail.push_back(cur);
        }

        int nc = g[cur][0];
        while (nc != cur) {
            cycle.push_back(nc);
            nc = g[nc][0];
        }
    }
}
