#include "gtest/gtest.h"
#include "../src/NRP/NRP.h"
#include "../src/NRP/NRPCycle.h"
#include "../src/NRP/NRPLine.h"
#include "../src/NRP/NRPtail.h"
#include <algorithm>

namespace nrp {
    class ContainNRPsTest : public ::testing::Test {
    protected:
        std::vector<aminoacid::Aminoacids::Aminoacid> amnacid;

        void genRandNrpPart(int partlen, nrpsprediction::NRPsPart& nrps_part, std::vector<std::vector<aminoacid::Aminoacids::Aminoacid>>& predict) {
            predict.resize(partlen);
            for (int j = 0; j < partlen; ++j) {
                std::vector<double> prob;
                std::vector<std::string> names;
                for (int g = 0; g < 3; ++g) {
                    prob.push_back(rand()%10 * 10);
                    predict[j].push_back(aminoacid::Aminoacids::Aminoacid(rand()%(int(aminoacid::Aminoacids::AMINOACID_CNT))));

                    names.push_back(aminoacid::Aminoacids::AMINOACID_NAMES[predict[j][predict[j].size() - 1]]);
                }

                std::sort(prob.rbegin(), prob.rend());
                std::stringstream ss;
                for (int g = 0; g < 3; ++g) {
                    ss << names[g] << "(" << prob[g] << ")";
                    if (g < 2) {
                        ss << ";";
                    }
                }

                for (int g = 0; g < 3; ++g) {
                    if (predict[j][g] == aminoacid::Aminoacids::Aminoacid::glu) {
                        predict[j].push_back(aminoacid::Aminoacids::me3_glu);
                    }

                    if (predict[j][g] == aminoacid::Aminoacids::Aminoacid::asn) {
                        predict[j].push_back(aminoacid::Aminoacids::OH_asn);
                    }
                }

                nrps_part.add_prediction(j + 1, ss.str());
            }
        }

        NRPCycle genRandCycleNRP(int len) {
            std::vector<std::string> strformula(len);
            std::vector<int> position(len);

            amnacid.resize(0);
            for (int i = 0; i < len; ++i) {
                position[i] = i;
                amnacid.push_back(aminoacid::Aminoacids::Aminoacid(rand() % aminoacid::Aminoacids::AMINOACID_CNT));
            }

            std::random_shuffle(position.begin(), position.end());

            NRPCycle res("", strformula, amnacid, position, "");
            res.aminoacids = amnacid;
            return res;
        }

        NRPLine genRandLineNRP(int len) {
            std::vector<std::string> strformula(len);
            std::vector<int> position(len);

            amnacid.resize(0);
            for (int i = 0; i < len; ++i) {
                position[i] = i;
                amnacid.push_back(aminoacid::Aminoacids::Aminoacid(rand() % aminoacid::Aminoacids::AMINOACID_CNT));
            }

            std::random_shuffle(position.begin(), position.end());

            NRPLine res("", strformula, amnacid, position, "");
            res.aminoacids = amnacid;
            return res;
        }

        bool eqval(int st, int delta, std::vector<std::vector<aminoacid::Aminoacids::Aminoacid>> predict) {
            for (int i = st, j = 0; j < predict.size(); ++j, i = (i + delta + amnacid.size()) % amnacid.size()) {
                bool has = false;

                for (int g = 0; g < predict[j].size(); ++g) {
                    if (aminoacid::Aminoacids::FORMULS[predict[j][g]] == aminoacid::Aminoacids::FORMULS[amnacid[i]]) {
                        has = true;
                    }
                }

                if (has == false) return false;
            }

            return true;
        }
    };


    //check if has segment in NRP then find it.
    TEST_F(ContainNRPsTest, findSegRandCycleTest) {
        for (int tst = 0; tst < 1000; ++tst) {
            int len = rand()%100 + 2;
            NRPCycle nrp = genRandCycleNRP(len);
            int bg = rand()%len, sz = rand()%len + 1;

            nrpsprediction::NRPsPart nrps_part("filename", "orf");

            int delta = 1;
            if (rand() % 2 == 0) {
                delta = -1;
            }

            int ed = (bg + (sz - 1) * delta + len)%len;

            int j = 0;
            for (int i = bg; j < sz; i = (i + delta + len)%len, ++j) {
                std::vector<double> prob;
                std::vector<std::string> names;
                for (int g = 0; g < 3; ++g) {
                    prob.push_back(rand()%10 * 10);
                    names.push_back(aminoacid::Aminoacids::AMINOACID_NAMES[rand()%(int(aminoacid::Aminoacids::AMINOACID_CNT))]);
                }

                names[rand()%3] = aminoacid::Aminoacids::AMINOACID_NAMES[int(amnacid[i])];
                std::sort(prob.rbegin(), prob.rend());
                std::stringstream ss;
                for (int g = 0; g < 3; ++g) {
                    ss << names[g] << "(" << prob[g] << ")";
                    if (g < 2) {
                        ss << ";";
                    }
                }
                nrps_part.add_prediction(j + 1, ss.str());
            }

            std::vector<NRP::Segment> segments = nrp.containNRPsPart(nrps_part);
            int rbg = bg;
            int red = ed;
            int rrev = 0;
            if (delta == 1) {
                if (red < rbg) {
                    red += len;
                }
            } else {
                rrev = 1;
                std::swap(rbg, red);
                if (red < rbg) {
                    red += len;
                }
            }


            bool hasRightAns = false;
            for (int i = 0; i < segments.size(); ++i) {
                if (segments[i].l == rbg && segments[i].r == red && segments[i].rev == rrev) {
                    hasRightAns = true;
                }
            }

            ASSERT_TRUE(hasRightAns);
        }
    }


    //check all find segment is right
    TEST_F(ContainNRPsTest, foundSegIsCorrectRandCycleTest) {
        for (int tst = 0; tst < 1000; ++tst) {
            int len = rand()%20 + 1;
            NRPCycle nrp = genRandCycleNRP(len);
            int partlen = rand()%len + 1;

            nrpsprediction::NRPsPart nrps_part("filename", "orf");
            std::vector<std::vector<aminoacid::Aminoacids::Aminoacid>> predict(partlen);

            genRandNrpPart(partlen, nrps_part, predict);

            std::vector<NRP::Segment> segments = nrp.containNRPsPart(nrps_part);
            for (int i = 0; i < segments.size(); ++i) {
                if (segments[i].rev == false) {
                    ASSERT_TRUE(eqval(segments[i].l, 1, predict));
                } else {
                    int r = segments[i].r;
                    if (segments[i].r >= len) {
                        r -= len;
                    }
                    ASSERT_TRUE(eqval(r, -1, predict));
                }
            }
        }
    }

    //check if has segment in NRP then find it.
    TEST_F(ContainNRPsTest, findSegRandLineTest) {
        for (int tst = 0; tst < 1000; ++tst) {
            int len = rand()%100 + 2;
            NRPLine nrp = genRandLineNRP(len);
            nrpsprediction::NRPsPart nrps_part("filename", "orf");
            int bg = rand()%len, ed = rand()%len;
            if (ed < bg) {
                std::swap(bg, ed);
            }
            int sz = ed - bg + 1;

            int delta = 1;
            if (rand() % 2 == 0) {
                delta = -1;
                std::swap(bg, ed);
            }

            int j = 0;
            for (int i = bg; j < sz; i += delta, ++j) {
                std::vector<double> prob;
                std::vector<std::string> names;
                for (int g = 0; g < 3; ++g) {
                    prob.push_back(rand()%10 * 10);
                    names.push_back(aminoacid::Aminoacids::AMINOACID_NAMES[rand()%(int(aminoacid::Aminoacids::AMINOACID_CNT))]);
                }

                names[rand()%3] = aminoacid::Aminoacids::AMINOACID_NAMES[int(amnacid[i])];
                std::sort(prob.rbegin(), prob.rend());
                std::stringstream ss;
                for (int g = 0; g < 3; ++g) {
                    ss << names[g] << "(" << prob[g] << ")";
                    if (g < 2) {
                        ss << ";";
                    }
                }
                nrps_part.add_prediction(j + 1, ss.str());
            }

            std::vector<NRP::Segment> segments = nrp.containNRPsPart(nrps_part);
            int rbg = std::min(bg, ed);
            int red = std::max(bg, ed);
            int rrev = (delta == -1);

            bool hasRightAns = false;
            for (int i = 0; i < segments.size(); ++i) {
                if (segments[i].l == rbg && segments[i].r == red && segments[i].rev == rrev) {
                    hasRightAns = true;
                }
            }

            ASSERT_TRUE(hasRightAns);
        }
    }


    //check all find segment is right
    TEST_F(ContainNRPsTest, foundSegIsCorrectRandLineTest) {
        for (int tst = 0; tst < 1000; ++tst) {
            int len = rand()%20 + 1;
            NRPLine nrp = genRandLineNRP(len);

            int partlen = rand()%len + 1;
            nrpsprediction::NRPsPart nrps_part ("filename", "orf");
            std::vector<std::vector<aminoacid::Aminoacids::Aminoacid>> predict(partlen);

            genRandNrpPart(partlen, nrps_part, predict);

            std::vector<NRP::Segment> segments = nrp.containNRPsPart(nrps_part);
            for (int i = 0; i < segments.size(); ++i) {
                ASSERT_TRUE(segments[i].r < len);
                ASSERT_TRUE(segments[i].l <= segments[i].r);
                if (segments[i].rev == false) {
                    ASSERT_TRUE(eqval(segments[i].l, 1, predict));
                } else {
                    ASSERT_TRUE(eqval(segments[i].r, -1, predict));
                }
            }
        }
    }
}

int main(int argc, char *argv[]) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}