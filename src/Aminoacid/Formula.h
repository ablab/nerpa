//
// Created by olga on 01.02.19.
//

#ifndef NRPSMATCHER_FORMULA_H
#define NRPSMATCHER_FORMULA_H

#include<string>
#include <map>
#include <cassert>

namespace aminoacid {
    class Formula {
    private:
        std::map<std::string, int> formula = {{"C", 0},
                                              {"H", 0},
                                              {"N", 0},
                                              {"O", 0},
                                              {"S", 0},
                                              {"Cl", 0}};
    public:
        Formula() = default;
        explicit Formula(std::string s) {
            std::string elem;
            int cnt_elem = 0;
            int sign = 1;
            for (char c : s) {
                if (c == '-') {
                    sign = -1;
                }
                if (c >= 'A' && c <= 'Z') {
                    if (cnt_elem == 0) {
                        cnt_elem = 1;
                    }
                    formula[elem] = sign * cnt_elem;
                    elem = c;
                    cnt_elem = 0;
                    sign = 1;
                } else if (c >= 'a' && c <= 'z') {
                    elem += c;
                } else if (c >= '0' && c <= '9') {
                    cnt_elem = cnt_elem * 10 + (c - '0');
                } else {
                    static_assert(true);
                }
            }

            if (cnt_elem == 0) {
                cnt_elem = 1;
            }
            formula[elem] = sign * cnt_elem;
        }

        bool operator == (const Formula& b) const {
            return formula == b.formula;
        }
    };
}


#endif //NRPSMATCHER_FORMULA_H
