#include "gtest/gtest.h"
#include "../src/NRP/NRP.h"
#include <algorithm>

namespace nrp {
    class ContainNRPsTest : public ::testing::Test {
    protected:
        std::vector<aminoacid::Aminoacids::Aminoacid> amnacid;

        NRP genRandNRP(int len) {
            amnacid.resize(0);
            for (int i = 0; i < len; ++i) {
                amnacid.push_back(aminoacid::Aminoacids::Aminoacid(rand() % aminoacid::Aminoacids::AMINOACID_CNT));
            }

            NRP res;
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


    TEST_F(ContainNRPsTest, trueRandTest) {
        for (int tst = 0; tst < 1000; ++tst) {
            int len = rand()%100 + 1;
            NRP nrp = genRandNRP(len);
            int bg = rand()%len, ed = rand()%len;

            nrpsprediction::NRPsPart nrps_part("filename", "orf");

            int delta = 1;
            if (rand() % 2 == 0) {
                delta = -1;
            }

            int j = 0;
            for (int i = bg; i != ed; i = (i + delta + len)%len, ++j) {
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

            ASSERT_TRUE(nrp.containNRPsPart(nrps_part));
        }
    }


    TEST_F(ContainNRPsTest, randTest) {
        for (int tst = 0; tst < 1000; ++tst) {
            int len = rand()%20 + 1;
            NRP nrp = genRandNRP(len);
            int partlen = rand()%len;

            nrpsprediction::NRPsPart nrps_part("filename", "orf");
            std::vector<std::vector<aminoacid::Aminoacids::Aminoacid>> predict(partlen);

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

                nrps_part.add_prediction(j + 1, ss.str());
            }

            bool ans = false;
            for (int l = 0; l < len; ++l) {
                if (eqval(l, 1, predict) || eqval(l, -1, predict)) {
                    ans = true;
                }
            }

            ASSERT_EQ(nrp.containNRPsPart(nrps_part), ans);
        }
    }


    TEST_F(ContainNRPsTest, minTest) {
        NRP nrp = genRandNRP(0);
        NRP nrp1 = genRandNRP(1);
        nrpsprediction::NRPsPart nrps_part("filename", "orf");

        ASSERT_FALSE(nrp.containNRPsPart(nrps_part)); //not sure what behavior should be expect
        ASSERT_TRUE(nrp1.containNRPsPart(nrps_part));
        nrps_part.add_prediction(1, "asn(90.0);asp(90.0);cys(50.0);leu(50.0);glu(50.0);gln(50.0);val(50.0);pip(50.0);gly(50.0);phe(50.0);ile(50.0);ala(50.0);hty(50.0);met(50.0);cha(50.0);orn(50.0)");

        ASSERT_FALSE(nrp.containNRPsPart(nrps_part));
    }


    TEST_F(ContainNRPsTest, maxTest) {
        //TODO
    }
}

int main(int argc, char *argv[]) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}