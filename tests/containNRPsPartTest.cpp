#include "gtest/gtest.h"
#include "../src/NRP/NRP.h"
#include "../src/NRP/NRPCycle.h"
#include "../src/NRP/NRPLine.h"
#include "../src/NRP/NRPtail.h"
#include <algorithm>
#include <cmath>
#include <Matcher/Matcher.h>
#include <Logger/log_writers.hpp>

namespace nrp {
    const double EPS = 1e-4;

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

            NRPCycle res("", strformula, amnacid, position, "", "");
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

            NRPLine res("", strformula, amnacid, position, "", "");
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

        void check_correct_score(matcher::Matcher::Match match, std::vector< nrpsprediction::NRPsPart> nrpParts, NRP::NRPType type) {
            std::vector<std::pair<int, int> > matchs = match.getMatchs();
            ASSERT_EQ(matchs.size(), amnacid.size());
            double score = 0;
            std::set<int> difseg;
            for (int i = 0 ; i < matchs.size(); ++i) {
                if (matchs[i].first == -1 && ((i == 0 && type != NRP::cycle) ||
                        (matchs[(i-1 + matchs.size())%matchs.size()].first != -1))) {
                    score -= 1;
                }

                if (matchs[i].first != -1) {
                    nrpsprediction::AminoacidPrediction amn_pred = nrpParts[matchs[i].first].getAminoacidsPrediction()[matchs[i].second];
                    score += amn_pred.getScore(amnacid[i]);
                    difseg.insert(matchs[i].first);
                }
            }
            score -= difseg.size();
            if (difseg.size() == 0 && type == NRP::cycle) {
                score -= 1;
            }
            ASSERT_LE(abs(score - match.score()), EPS);
        }

        void check_correct_seq(matcher::Matcher::Match match, std::vector< nrpsprediction::NRPsPart> nrpParts) {
            std::vector<std::vector<int> > post(nrpParts.size());
            for (int i = 0; i < post.size(); ++i) {
                post[i].resize(nrpParts[i].getAminoacidsPrediction().size(), -1);
            }

            std::vector<std::pair <int, int> > matchs = match.getMatchs();
            for (int i = 0; i < matchs.size(); ++i) {
                if (matchs[i].first != -1) {
                    ASSERT_LE(matchs[i].first, nrpParts.size() - 1);
                    ASSERT_GE(matchs[i].second, 0);
                    ASSERT_LE(matchs[i].second, post[matchs[i].first].size() - 1);
                    post[matchs[i].first][matchs[i].second] = i;
                }
            }

            for (int i = 0; i < post.size(); ++i) {
                bool isOne = true, isUp = true, isDown = true;
                for (int j = 0; j < post[i].size(); ++j) {
                    if (post[i][j] == -1) {
                        isUp = false;
                        isDown = false;
                    }
                    if (post[i][j] != -1) {
                        isOne = false;
                    }
                    if (j != 0 && post[i][j] + 1 != post[i][j - 1]) {
                        isUp = false;
                    }
                    if (j != 0 && post[i][j] != post[i][j - 1] + 1) {
                        isDown = false;
                    }
                }
                ASSERT_TRUE(isOne || isUp || isDown);
            }
        }

        void check_correct_cycle_seq(matcher::Matcher::Match match, std::vector< nrpsprediction::NRPsPart> nrpParts) {
            std::vector<std::vector<int> > post(nrpParts.size());
            for (int i = 0; i < post.size(); ++i) {
                post[i].resize(nrpParts[i].getAminoacidsPrediction().size(), -1);
            }

            std::vector<std::pair <int, int> > matchs = match.getMatchs();
            for (int i = 0; i < matchs.size(); ++i) {
                if (matchs[i].first != -1) {
                    ASSERT_LE(matchs[i].first, nrpParts.size() - 1);
                    ASSERT_GE(matchs[i].second, 0);
                    ASSERT_LE(matchs[i].second, post[matchs[i].first].size() - 1);
                    post[matchs[i].first][matchs[i].second] = i;
                }
            }

            for (int i = 0; i < post.size(); ++i) {
                bool isOne = true, isUp = true, isDown = true;
                for (int j = 0; j < post[i].size(); ++j) {
                    if (post[i][j] == -1) {
                        isUp = false;
                        isDown = false;
                    }
                    if (post[i][j] != -1) {
                        isOne = false;
                    }
                    if (j != 0 && (post[i][j] + 1)%amnacid.size() != post[i][j - 1]) {
                        isUp = false;
                    }
                    if (j != 0 && post[i][j] != (post[i][j - 1] + 1)%amnacid.size()) {
                        isDown = false;
                    }
                }
                ASSERT_TRUE(isOne || isUp || isDown);
            }
        }


        void check_seq_match_with_nrp(matcher::Matcher::Match match, std::vector< nrpsprediction::NRPsPart> nrpParts) {
            std::vector<std::pair <int, int> > matchs = match.getMatchs();
            for (int i = 0; i < matchs.size(); ++i) {
                if (matchs[i].first != -1) {
                    nrpsprediction::AminoacidPrediction apr = nrpParts[matchs[i].first].getAminoacidsPrediction()[matchs[i].second];
                    ASSERT_TRUE(apr.contain(amnacid[i]));
                }
            }
        }

        nrpsprediction::NRPsPart getSubPart(int bg, int sz, int delta, double &score) {
            nrpsprediction::NRPsPart nrps_part("filename", "orf");
            int len = amnacid.size();
            int j = 0;
            for (int i = bg; j < sz; i = (i + delta + len)%len, ++j) {
                std::vector<double> prob;
                std::vector<std::string> names;
                for (int g = 0; g < 3; ++g) {
                    prob.push_back(60 + rand()%4 * 10);
                    names.push_back(aminoacid::Aminoacids::AMINOACID_NAMES[rand()%(int(aminoacid::Aminoacids::AMINOACID_CNT))]);
                }

                int right_AA_pos = rand()%3;
                names[right_AA_pos] = aminoacid::Aminoacids::AMINOACID_NAMES[int(amnacid[i])];
                std::sort(prob.rbegin(), prob.rend());
                std::stringstream ss;
                for (int g = 0; g < 3; ++g) {
                    ss << names[g] << "(" << prob[g] << ")";
                    if (g < 2) {
                        ss << ";";
                    }
                }
                nrps_part.add_prediction(j + 1, ss.str());
                auto predictions = nrps_part.getAminoacidsPrediction();
                score += predictions.back().getScore(amnacid[i]);
            }

            return nrps_part;
        }
    };


    //check if has segment in NRP then find it.
    TEST_F(ContainNRPsTest, containTrueRandCycleTest) {
        double scr = 0;
        for (int tst = 0; tst < 1000; ++tst) {
            int len = rand()%100 + 2;
            NRPCycle nrp = genRandCycleNRP(len);
            int bg = rand()%len, sz = rand()%len + 1;

            int delta = 1;
            if (rand() % 2 == 0) {
                delta = -1;
            }

            int ed = (bg + (sz - 1) * delta + len)%len;

            nrpsprediction::NRPsPart nrps_part = getSubPart(bg, sz, delta, scr);

            std::vector<nrpsprediction::NRPsPart> parts;
            parts.push_back(nrps_part);
            matcher::Matcher matcher1(nrp, nrpsprediction::NRPsPrediction(parts));

            std::vector<NRP::Segment> segments = matcher1.matche_seg(nrps_part);
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
    TEST_F(ContainNRPsTest, containRandCycleTest) {
        for (int tst = 0; tst < 1000; ++tst) {
            int len = rand()%20 + 1;
            NRPCycle nrp = genRandCycleNRP(len);
            int partlen = rand()%len + 1;

            nrpsprediction::NRPsPart nrps_part("filename", "orf");
            std::vector<std::vector<aminoacid::Aminoacids::Aminoacid>> predict(partlen);

            genRandNrpPart(partlen, nrps_part, predict);

            std::vector<nrpsprediction::NRPsPart> parts;
            parts.push_back(nrps_part);
            matcher::Matcher matcher1(nrp, nrpsprediction::NRPsPrediction(parts));

            std::vector<NRP::Segment> segments = matcher1.matche_seg(nrps_part);
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
    TEST_F(ContainNRPsTest, containTrueRandLineTest) {
        double scr = 0;
        for (int tst = 0; tst < 1000; ++tst) {
            int len = rand()%100 + 2;
            NRPLine nrp = genRandLineNRP(len);

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

            nrpsprediction::NRPsPart nrps_part = getSubPart(bg, sz, delta, scr);

            std::vector<nrpsprediction::NRPsPart> parts;
            parts.push_back(nrps_part);
            matcher::Matcher matcher1(nrp, nrpsprediction::NRPsPrediction(parts));

            std::vector<NRP::Segment> segments = matcher1.matche_seg(nrps_part);
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
    TEST_F(ContainNRPsTest, containRandLineTest) {
        for (int tst = 0; tst < 1000; ++tst) {
            int len = rand()%20 + 1;
            NRPLine nrp = genRandLineNRP(len);

            int partlen = rand()%len + 1;
            nrpsprediction::NRPsPart nrps_part ("filename", "orf");
            std::vector<std::vector<aminoacid::Aminoacids::Aminoacid>> predict(partlen);

            genRandNrpPart(partlen, nrps_part, predict);

            std::vector<nrpsprediction::NRPsPart> parts;
            parts.push_back(nrps_part);
            matcher::Matcher matcher1(nrp, nrpsprediction::NRPsPrediction(parts));

            std::vector<NRP::Segment> segments = matcher1.matche_seg(nrps_part);
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

    TEST_F(ContainNRPsTest, coverRandLineTest) {
        for (int tst = 0; tst < 1000; ++tst) {
            int len = rand()%20 + 1;
            NRPLine nrp = genRandLineNRP(len);

            int cntbp = rand()%5 + 1;
            std::vector<int> bps(cntbp);
            for (int i = 0; i < cntbp; ++i) {
                bps[i] = rand()%len;
            }

            bps.push_back(0);
            bps.push_back(len);

            std::sort(bps.begin(), bps.end());

            std::vector<nrpsprediction::NRPsPart> nrpParts;
            double res_score = 0;
            for (int i = 1; i < bps.size(); ++i) {
                if (bps[i] != bps[i - 1]) {
                    if (rand() % 2 == 0) {
                        if (rand() % 2 == 0) {
                            nrpParts.push_back(getSubPart(bps[i - 1], (bps[i] - bps[i - 1]), 1, res_score));
                        } else {
                            nrpParts.push_back(getSubPart(bps[i] - 1, (bps[i] - bps[i - 1]), -1, res_score));
                        }
                    }
                    res_score -= 1;
                }
            }

            for (int i = 0; i < 10; ++i) {
                int partlen = rand()%10 + 1;
                nrpsprediction::NRPsPart nrps_part ("filename", "orf");
                std::vector<std::vector<aminoacid::Aminoacids::Aminoacid>> predict(partlen);

                genRandNrpPart(partlen, nrps_part, predict);
                nrpParts.push_back(nrps_part);
            }

            std::random_shuffle(nrpParts.begin(), nrpParts.end());
            nrpsprediction::NRPsPrediction nrpsPrediction(nrpParts);

            matcher::Matcher matcher(nrp, nrpsPrediction);
            matcher::Matcher::Match match = matcher.getMatch();

            ASSERT_GE(match.score() - res_score, -EPS);
            check_correct_score(match, nrpParts, NRP::line);
            check_correct_seq(match, nrpParts);
            check_seq_match_with_nrp(match, nrpParts);
        }
    }

    TEST_F(ContainNRPsTest, coverRandCycleTest) {
        for (int tst = 0; tst < 1000; ++tst) {
            int len = rand()%20 + 1;
            NRPCycle nrp = genRandCycleNRP(len);

            int cntbp = rand()%5 + 1;
            std::vector<int> bps(cntbp);
            for (int i = 0; i < cntbp; ++i) {
                bps[i] = rand()%len;
            }

            std::sort(bps.begin(), bps.end());

            std::vector<nrpsprediction::NRPsPart> nrpParts;
            double res_score = 0;
            for (int i = 0; i < bps.size(); ++i) {
                int pi = (i - 1 + bps.size()) % bps.size();
                if (bps[i] != bps[pi] || (bps[i] == bps[pi] && i == 0)) {
                    res_score -= 1;
                    if (rand() % 2 == 0) {
                        if (rand() % 2 == 0) {
                            nrpParts.push_back(getSubPart(bps[pi], (bps[i] - bps[pi] + len) % len, 1, res_score));
                        } else {
                            nrpParts.push_back(getSubPart((bps[i] - 1 + len) % len, (bps[i] - bps[pi] + len)%len, -1, res_score));
                        }
                    }
                }
            }

            for (int i = 0; i < 10; ++i) {
                int partlen = rand()%10 + 1;
                nrpsprediction::NRPsPart nrps_part ("filename", "orf");
                std::vector<std::vector<aminoacid::Aminoacids::Aminoacid>> predict(partlen);

                genRandNrpPart(partlen, nrps_part, predict);
                nrpParts.push_back(nrps_part);
            }

            std::random_shuffle(nrpParts.begin(), nrpParts.end());
            nrpsprediction::NRPsPrediction nrpsPrediction(nrpParts);

            matcher::Matcher matcher(nrp, nrpsPrediction);
            matcher::Matcher::Match match = matcher.getMatch();

            ASSERT_GE(match.score() - res_score, -EPS);
            check_correct_score(match, nrpParts, NRP::cycle);
            check_correct_cycle_seq(match, nrpParts);
            check_seq_match_with_nrp(match, nrpParts);
        }
    }
}

int main(int argc, char *argv[]) {
    ::testing::InitGoogleTest(&argc, argv);
    logging::create_logger("");
    return RUN_ALL_TESTS();
}