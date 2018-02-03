#ifndef NRPSMATCHER_NRPBUILDER_H
#define NRPSMATCHER_NRPBUILDER_H

#include <string>
#include "NRP.h"

namespace nrp {
    class NRPBuilder {
    private:
        static std::vector<aminoacid::Aminoacids::Aminoacid> aminoacids_by_pos(
                const std::vector<aminoacid::Aminoacids::Aminoacid> &aminoacids, const std::vector<int> &pos);
    public:
        static NRP* build(std::string fragment_graph);

        static bool isCycle(std::vector<std::vector<int>> &g, std::vector<std::vector<int>> &gr);

        static bool isLine(std::vector<std::vector<int>> &g, std::vector<std::vector<int>> &gr);

        static bool isTail(std::vector<std::vector<int>> &g, std::vector<std::vector<int>> &gr);

        static std::vector<int> parseCycle(std::vector<std::vector<int>> &g, std::vector<std::vector<int>> &gr);

        static std::vector<int> parseLine(std::vector<std::vector<int>> &g, std::vector<std::vector<int>> &gr);

        static void parseTail(std::vector<std::vector<int>> &g, std::vector<std::vector<int>> &gr, std::vector<int> &tail,
                              std::vector<int> &cycle);
    };
}

#endif //NRPSMATCHER_NRPBUILDER_H
