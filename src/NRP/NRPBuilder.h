#ifndef NRPSMATCHER_NRPBUILDER_H
#define NRPSMATCHER_NRPBUILDER_H

#include <string>
#include "NRP.h"

namespace nrp {
    class NRPBuilder {
    private:
        enum Elem{C, H, Cl, N, O, S, ELEM_CNT};
        static const std::string ELEM_NAME[ELEM_CNT];

        static std::vector<aminoacid::Aminoacids::Aminoacid> aminoacids_by_pos(
                const std::vector<aminoacid::Aminoacids::Aminoacid> &aminoacids, const std::vector<int> &pos);

        static int get_elem_id(std::string c);
        static std::vector<int> parse_formula(std::string formula);
        static std::string to_string(std::vector<int> formula);

        static bool isCycle(std::vector<std::vector<int>> &g, std::vector<std::vector<int>> &gr);

        static bool isLine(std::vector<std::vector<int>> &g, std::vector<std::vector<int>> &gr);

        static bool isTail(std::vector<std::vector<int>> &g, std::vector<std::vector<int>> &gr);

        static std::vector<int> parseCycle(std::vector<std::vector<int>> &g, std::vector<std::vector<int>> &gr);

        static std::vector<int> parseLine(std::vector<std::vector<int>> &g, std::vector<std::vector<int>> &gr);

        static void parseTail(std::vector<std::vector<int>> &g, std::vector<std::vector<int>> &gr, std::vector<int> &tail,
                              std::vector<int> &cycle);

    public:
        static NRP* build(std::string fragment_graph, std::string extra_info);

        static bool isConnected(std::vector<std::vector<int>> &g, std::vector<std::vector<int>> &gr);

        static void dfs(int v, std::vector<int> &used, std::vector<std::vector<int>> &g, std::vector<std::vector<int>> &gr);
    };
}

#endif //NRPSMATCHER_NRPBUILDER_H
