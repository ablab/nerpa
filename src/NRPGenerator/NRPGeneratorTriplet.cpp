#include <NRP/NRPCycle.h>
#include <NRP/NRPLine.h>
#include <algorithm>
#include <NRP/NRPtail.h>
#include "NRPGeneratorTriplet.h"

namespace nrp_generator {
    NRPGeneratorTriplet::NRPGeneratorTriplet(const std::vector<nrp::NRP *> &mols) : NRPGenerator(mols) {
        graph.resize(aminoacid::Aminoacids::AMINOACID_CNT*aminoacid::Aminoacids::AMINOACID_CNT,
                     std::vector<int> (aminoacid::Aminoacids::AMINOACID_CNT*aminoacid::Aminoacids::AMINOACID_CNT, 0));

        cntVert.resize(aminoacid::Aminoacids::AMINOACID_CNT*aminoacid::Aminoacids::AMINOACID_CNT, 0);
        for (auto mol : mols) {
            int mollen = mol->getLen();
            if (mollen < 2) continue;
            for (int j = 2; j < mollen + 2; ++j) {
                graph[getVertexId(mol->getAminoacid(j - 2), mol->getAminoacid((j - 1)%mollen))]
                [getVertexId(mol->getAminoacid((j - 1)%mollen), mol->getAminoacid(j%mollen))] += 1;
                cntVert[getVertexId(mol->getAminoacid((j - 1)%mollen), mol->getAminoacid(j%mollen))] += 1;
            }
        }


        for (int i = 0; i < graph.size(); ++i) {
            for (int j = 1; j < graph[i].size(); ++j) {
                graph[i][j] += graph[i][j - 1];
            }
        }
        for (int i = 1; i < cntVert.size(); ++i) {
            cntVert[i] += cntVert[i - 1];
        }
    }

    nrp::NRP *NRPGeneratorTriplet::generate(nrp::NRP::NRPType type, int len) {
        if (type == nrp::NRP::branch_cycle) {
            return generateNRPBranchCycle(len);
        }
        if (type == nrp::NRP::cycle) {
            return generateNRPCycle(len);
        }
        if (type == nrp::NRP::line) {
            return generateNRPLine(len);
        }
    }

    nrp::NRP *NRPGeneratorTriplet::generateNRPCycle(int len) {
        std::vector<std::string> strformula(len);
        std::vector<int> position(len, 0);
        std::vector<aminoacid::Aminoacids::Aminoacid> amnacid;

        amnacid.resize(0);
        int strv = genStartVert();
        amnacid.push_back(getVertexAA(strv).first);
        amnacid.push_back(getVertexAA(strv).second);
        position[1] = 1;
        for (int i = 2; i < len; ++i) {
            position[i] = i;
            amnacid.push_back(genAA(strv));
        }

        nrp::NRP *res = new nrp::NRPCycle("", strformula, amnacid, position, "", "");
        return res;
    }

    nrp::NRP *NRPGeneratorTriplet::generateNRPLine(int len) {
        std::vector<aminoacid::Aminoacids::Aminoacid> amnacid;
        std::vector<std::string> strformula(len);
        std::vector<int> position(len);

        amnacid.resize(0);
        int strv = genStartVert();
        amnacid.push_back(getVertexAA(strv).first);
        amnacid.push_back(getVertexAA(strv).second);
        position[1] = 1;
        for (int i = 2; i < len; ++i) {
            position[i] = i;
            amnacid.push_back(genAA(strv));
        }

        nrp::NRP *res = new nrp::NRPLine("", strformula, amnacid, position, "", "");
        return res;
    }

    nrp::NRP *NRPGeneratorTriplet::generateNRPBranchCycle(int len) {
        int tailLen = rand() % len;
        std::vector<std::string> strformula(len);
        std::vector<int> position(len);
        std::vector<aminoacid::Aminoacids::Aminoacid> amnacid;

        amnacid.resize(0);
        int strv = genStartVert();
        amnacid.push_back(getVertexAA(strv).first);
        amnacid.push_back(getVertexAA(strv).second);
        position[1] = 1;
        for (int i = 2; i < len; ++i) {
            position[i] = i;
            amnacid.push_back(genAA(strv));
        }
        nrp::NRPLine res1 = nrp::NRPLine("", strformula, amnacid, position, "", "");

        std::reverse(position.begin() + tailLen, position.end());
        nrp::NRPLine res2 = nrp::NRPLine("", strformula, amnacid, position, "", "");
        nrp::NRP *res = new nrp::NRPtail(res1, res2);
        return res;
    }

    Aminoacid NRPGeneratorTriplet::genAA(int& v) {
        int rv = rand()%graph[v][graph[v].size() - 1];
        int u = 0;
        while (graph[v][u] < rv) {
            ++u;
        }
        v = u;
        return getVertexAA(v).second;
    }

    int NRPGeneratorTriplet::genStartVert() {
        int rv = rand()%cntVert[cntVert.size() - 1];
        int v = 0;
        while (cntVert[v] < rv) {
            ++v;
        }
        return v;
    }

    int NRPGeneratorTriplet::getVertexId(Aminoacid aa1, Aminoacid aa2) {
        return (int)aa1 * aminoacid::Aminoacids::AMINOACID_CNT + (int)aa2;
    }

    std::pair<Aminoacid, Aminoacid> NRPGeneratorTriplet::getVertexAA(int id) {
        return std::make_pair(Aminoacid(id/aminoacid::Aminoacids::AMINOACID_CNT),
                              Aminoacid(id%aminoacid::Aminoacids::AMINOACID_CNT));
    }
}