#include <algorithm>
#include "NRPGenerator.h"
#include "NRP/NRPCycle.h"
#include "NRP/NRPLine.h"
#include "NRP/NRPtail.h"

namespace nrp_generator {
    using namespace nrp;
    NRPGenerator::NRPGenerator(std::vector<nrp::NRP *> mols) {
        cntAA.resize(aminoacid::Aminoacid::AMINOACID_CNT, 0);
        for (int i = 0; i < mols.size(); ++i) {
            int mollen = mols[i]->getLen();
            for (int j = 0; j < mollen; ++j) {
                ++sumAA;
                cntAA[mols[i]->getAminoacid(j)] += 1;
            }
        }
    }

    nrp::NRP *NRPGenerator::generate(nrp::NRP::NRPType type, int len) {
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

    nrp::NRP *NRPGenerator::generateNRPCycle(int len) {
        std::vector<std::string> strformula(len);
        std::vector<int> position(len);
        std::vector<aminoacid::Aminoacid::AminoacidId> amnacid;

        amnacid.resize(0);
        for (int i = 0; i < len; ++i) {
            position[i] = i;
            amnacid.push_back(genAA());
        }

        nrp::NRP *res = new NRPCycle("", strformula, amnacid, position, "", "");
        return res;
    }

    nrp::NRP *NRPGenerator::generateNRPLine(int len) {
        std::vector<aminoacid::Aminoacid::AminoacidId> amnacid;
        std::vector<std::string> strformula(len);
        std::vector<int> position(len);

        amnacid.resize(0);
        for (int i = 0; i < len; ++i) {
            position[i] = i;
            amnacid.push_back(genAA());
        }

        nrp::NRP *res = new NRPLine("", strformula, amnacid, position, "", "");
        return res;
    }

    nrp::NRP *NRPGenerator::generateNRPBranchCycle(int len) {
        int tailLen = rand() % len;
        std::vector<std::string> strformula(len);
        std::vector<int> position(len);
        std::vector<aminoacid::Aminoacid::AminoacidId> amnacid;

        amnacid.resize(0);
        for (int i = 0; i < len; ++i) {
            position[i] = i;
            amnacid.push_back(genAA());
        }
        nrp::NRPLine res1 = NRPLine("", strformula, amnacid, position, "", "");

        std::reverse(position.begin() + tailLen, position.end());
        nrp::NRPLine res2 = NRPLine("", strformula, amnacid, position, "", "");
        nrp::NRP *res = new nrp::NRPtail(res1, res2);
        return res;
    }

    aminoacid::Aminoacid::AminoacidId NRPGenerator::genAA() {
        int rv = rand() % sumAA;
        int curSum = 0;
        for (int i = 0; i < aminoacid::Aminoacid::AMINOACID_CNT; ++i, curSum += cntAA[i]) {
            if (curSum <= rv && rv < curSum + cntAA[i]) {
                return aminoacid::Aminoacid::AminoacidId(i);
            }
        }

        return aminoacid::Aminoacid::AMINOACID_CNT;
    }
}