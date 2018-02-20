#include <fstream>
#include <sstream>
#include <iostream>
#include <cassert>
#include "NRPBuilder.h"
#include "NRPCycle.h"
#include "NRPLine.h"
#include "NRPtail.h"

const std::string nrp::NRPBuilder::ELEM_NAME[ELEM_CNT] = {"C", "H", "Cl", "N", "O", "S"};

nrp::NRP* nrp::NRPBuilder::build(std::string fragment_graph, std::string extra_info) {
    std::ifstream in(fragment_graph);
    std::string s;
    std::vector<std::string> strformula;
    std::vector<std::vector<int> > formuls;
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
        formuls.push_back(parse_formula(formula));
    }

    std::vector<std::vector<int> > g(formuls.size());
    std::vector<std::vector<int> > gr(formuls.size());

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
                formuls[b][int(Elem::H)] += 1;
                formuls[e][int(Elem::H)] += 1;
                formuls[e][int(Elem::O)] += 1;
            }
        }
    }

    for (int i = 0; i < formuls.size(); ++i) {
        aminoacids.push_back(aminoacid::Aminoacids::get_aminoacid_from_formula(to_string(formuls[i])));
    }

    in.close();

    if (!isConnected(g, gr)) {
        return nullptr;
    }

    if (isCycle(g, gr)) {
        std::vector<int> pos = parseCycle(g, gr);
        std::vector<aminoacid::Aminoacids::Aminoacid> resaacid = aminoacids_by_pos(aminoacids, pos);

        return new nrp::NRPCycle(fragment_graph, strformula, resaacid, pos, graph, extra_info);
    } else if (isLine(g, gr)) {
        std::vector<int> pos = parseLine(g, gr);
        std::vector<aminoacid::Aminoacids::Aminoacid> resaacid = aminoacids_by_pos(aminoacids, pos);

        return new NRPLine(fragment_graph, strformula, resaacid, pos, graph, extra_info);
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

        std::vector<aminoacid::Aminoacids::Aminoacid> resaacid1 = aminoacids_by_pos(aminoacids, pos1);
        std::vector<aminoacid::Aminoacids::Aminoacid> resaacid2 = aminoacids_by_pos(aminoacids, pos2);

        NRPLine ver1(fragment_graph, strformula, resaacid1, pos1, graph, extra_info), ver2(fragment_graph, strformula, resaacid2, pos2, graph, extra_info);

        return new NRPtail(ver1, ver2);
    }

    return nullptr;
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
        if (g[i].size() == 0 && gr[i].size() != 1) {
            return false;
        }
        if (g[i].size() == 0 && gr[i].size() == 1) {
            cnt01++;
        }
        if (g[i].size() > 1) {
            return false;
        }
    }

    for (int i = 0; i < (int)gr.size(); ++i) {
        if (gr[i].size() == 0 && g[i].size() != 1) {
            return false;
        }
        if (gr[i].size() == 0 && g[i].size() == 1) {
            cnt02++;
        }
        if (gr[i].size() > 1) {
            return false;
        }
    }

    return (cnt01 == 1 && cnt02 == 1);
}

bool nrp::NRPBuilder::isTail(std::vector<std::vector<int>> &g, std::vector<std::vector<int>> &gr) {
    int cnt01 = 0, cnt02 = 0;
    int cnt2 = 0;
    for (int i = 0; i < (int)g.size(); ++i) {
        if (g[i].size() == 0) {
            cnt01++;
        }
        if (g[i].size() > 2) {
            return false;
        }
        if (g[i].size() == 2) {
            cnt2 += 1;
        }
    }

    for (int i = 0; i < (int)gr.size(); ++i) {
        if (gr[i].size() == 0) {
            cnt02++;
        }
        if (gr[i].size() > 2) {
            return false;
        }
        if (gr[i].size() == 2) {
            cnt2 += 1;
        }
    }

    return (cnt01 + cnt02 == 1 && cnt2 == 1);
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
        while (g[cur].size() != 2) {
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
        while (gr[cur].size() != 2) {
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

int nrp::NRPBuilder::get_elem_id(std::string c) {
    for (int i = 0; i < ELEM_CNT; ++i) {
        if (ELEM_NAME[i] == c) {
            return i;
        }
    }

    return ELEM_CNT;
    assert(false);
}

std::string nrp::NRPBuilder::to_string(std::vector<int> formula) {
    std::stringstream res;
    for (int i = 0; i < formula.size(); ++i) {
        if (formula[i] != 0) {
            res << ELEM_NAME[i];
        }
        if (formula[i] > 1) {
            res << formula[i];
        }
    }

    return res.str();
}

std::vector<int> nrp::NRPBuilder::parse_formula(std::string formula) {
    std::vector<int> res(ELEM_CNT, 0);
    Elem curelem = ELEM_CNT;
    for (int i = 0; i < formula.size(); ++i) {
        if (formula[i] >= 'A' && formula[i] <= 'Z') {
            if (curelem != ELEM_CNT) {
                if (res[curelem] == 0) {
                    ++res[curelem];
                }
            }

            std::string nm = "";
            nm += formula[i];
            if (i + 1 < formula.size() && formula[i + 1] >= 'a' && formula[i + 1] <= 'z') {
                nm += formula[i + 1];
                ++i;
            }
            curelem = Elem(get_elem_id(nm));
            if (curelem == ELEM_CNT) {
                res[0] = 179;
                return res;
            }
        } else {
            res[curelem] = res[curelem] * 10 + (formula[i] - '0');
        }
    }

    if (curelem != ELEM_CNT) {
        if (res[curelem] == 0) {
            ++res[curelem];
        }
    }

    return res;
}

bool nrp::NRPBuilder::isConnected(std::vector<std::vector<int>> &g, std::vector<std::vector<int>> &gr) {
    std::vector<int> used(g.size(), 0);
    dfs(0, used, g, gr);
    for (int i = 0; i < g.size(); ++i) {
        if (used[i] == 0) {
            return false;
        }
    }

    return true;
}

void nrp::NRPBuilder::dfs(int v, std::vector<int> &used, std::vector<std::vector<int>> &g,
                          std::vector<std::vector<int>> &gr) {
    used[v] = 1;
    for (int i = 0; i < g[v].size(); ++i) {
        if (used[g[v][i]] == 0) {
            dfs(g[v][i], used, g, gr);
        }
    }


    for (int i = 0; i < gr[v].size(); ++i) {
        if (used[gr[v][i]] == 0) {
            dfs(gr[v][i], used, g, gr);
        }
    }
}
