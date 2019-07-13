//
// Created by olga on 01.02.19.
//

#ifndef NRPSMATCHER_FORMULA_H
#define NRPSMATCHER_FORMULA_H

#include <string>
#include <map>
#include <cassert>
#include <iostream>
#include <sstream>
#include <vector>

namespace aminoacid {
    class Formula {
    private:
        std::map<std::string, int> formula = {{"", 0},
                                              {"C", 0},
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
                    if (elem != "") formula[elem] = sign * cnt_elem;
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

            if (elem != "") formula[elem] = sign * cnt_elem;
        }

        bool operator == (const Formula& b) const {
            return formula == b.formula;
        }

        Formula operator - (const Formula& b) const {
            Formula res;
            res.formula = formula;
            for (auto const& x: b.formula) {
                if (res.formula.count(x.first)) {
                    res.formula[x.first] -= x.second;
                } else {
                    res.formula[x.first] = -x.second;
                }
            }

            return res;
        }

        Formula& operator += (const Formula& b) {
            for (auto const& x: b.formula) {
                if (formula.count(x.first)) {
                    formula[x.first] += x.second;
                } else {
                    formula[x.first] = x.second;
                }
            }
            return  *this;
        }

        void print() const {
            for (auto const& x: formula) {
                std::cerr << "formula[" << x.first << "] = " << x.second << "\n";
            }
        }

        std::string toString() {
            std::vector<std::string> elems({"C","H","N","O","S","Cl"});
            std::stringstream ss;
            for (auto elem : elems) {
                if (formula[elem] != 0) {
                    ss << elem;
                    if (formula[elem] != 1) {
                        ss << formula[elem];
                    }
                }
            }

            return ss.str();
        }
    };
}


#endif //NRPSMATCHER_FORMULA_H
