#include "gtest/gtest.h"
#include "../src/NRP/NRP.h"
#include "../src/NRP/NRPCycle.h"
#include "../src/NRP/NRPLine.h"
#include "../src/NRP/NRPBranch.h"
#include <algorithm>
#include <cmath>
#include <Matcher/Matcher.h>
#include <Logger/log_writers.hpp>
#include <boost/concept_check.hpp>
#include <NRPsPrediction/Builders/Nrpspredictor2Builder.h>
#include <Matcher/Score/Base/ScoreFullMatch.h>
#include <memory>

namespace nrp {
    typedef matcher::Segment Segment;
    const double EPS = 1e-4;

    class ContainNRPsTest : public ::testing::Test {
    protected:
        std::vector<aminoacid::Aminoacid> amnacid;
        virtual void SetUp() {
            aminoacid::AminoacidInfo::init("../../resources/aminoacids.tsv", "NRPSPREDICTOR2");
        }

        void genRandNrpPart(int partlen, nrpsprediction::OrfPrediction& nrps_part,
                            std::vector<std::vector<aminoacid::Aminoacid>>& predict) {
            predict.resize(partlen);
            for (int j = 0; j < partlen; ++j) {
                std::vector<double> prob;
                std::vector<std::string> names;
                for (int g = 0; g < 3; ++g) {
                    prob.push_back(rand()%10 * 10);
                    predict[j].push_back(aminoacid::Aminoacid(rand()%(aminoacid::AminoacidInfo::AMINOACID_CNT - 1)));

                    names.push_back(predict[j][predict[j].size() - 1].get_name());
                }

                std::sort(prob.rbegin(), prob.rend());
                std::stringstream ss;
                for (int g = 0; g < 3; ++g) {
                    ss << names[g] << "(" << prob[g] << ")";
                    if (g < 2) {
                        ss << ";";
                    }
                }

                /*for (int g = 0; g < 3; ++g) {
                    if (predict[j][g] == aminoacid::Aminoacid::AminoacidId::glu) {
                        predict[j].push_back(aminoacid::Aminoacid::me3_glu);
                    }

                    if (predict[j][g] == aminoacid::Aminoacid::AminoacidId::asn) {
                        predict[j].push_back(aminoacid::Aminoacid::OH_asn);
                    }
                }*/

                using namespace nrpsprediction;
                nrps_part.add_prediction(j + 1, AAdomainPrediction(j + 1,
                                                                   Nrpspredictor2Builder::parse_predictions(ss.str())));
            }
        }

        std::shared_ptr<nrp::NRP> genRandCycleNRP(int len) {
            std::vector<std::string> strformula(len);
            std::vector<int> position(len);

            amnacid.resize(0);
            for (int i = 0; i < len; ++i) {
                position[i] = i;
                amnacid.push_back(aminoacid::Aminoacid(rand()%(aminoacid::AminoacidInfo::AMINOACID_CNT - 1)));
            }

            std::random_shuffle(position.begin(), position.end());

            std::shared_ptr<nrp::NRP> res = std::make_shared<nrp::NRPCycle>("", strformula, amnacid, position, "", "");
            res->aminoacids = amnacid;
            return res;
        }

        std::shared_ptr<nrp::NRP> genRandLineNRP(int len) {
            std::vector<std::string> strformula(len);
            std::vector<int> position(len);

            amnacid.resize(0);
            for (int i = 0; i < len; ++i) {
                position[i] = i;
                amnacid.push_back(aminoacid::Aminoacid(rand()%(aminoacid::AminoacidInfo::AMINOACID_CNT - 1)));
            }

            std::random_shuffle(position.begin(), position.end());

            std::shared_ptr<nrp::NRP> res = std::make_shared<nrp::NRPLine>("", strformula, amnacid, position, "", "");
            res->aminoacids = amnacid;
            return res;
        }

        bool eqval(int st, int delta, std::vector<std::vector<aminoacid::Aminoacid>> predict) {
            for (int i = st, j = 0; j < predict.size(); ++j, i = (i + delta + amnacid.size()) % amnacid.size()) {
                bool has = false;

                for (int g = 0; g < predict[j].size(); ++g) {
                    if (predict[j][g] == amnacid[i]) {
                        has = true;
                    }
                }

                if (has == false) return false;
            }

            return true;
        }

        void check_correct_score(matcher::Matcher::Match match, std::vector< nrpsprediction::OrfPrediction> nrpParts, NRP::NRPType type) {
            nrpsprediction::BgcPrediction prediction(nrpParts);
            std::vector<std::pair<int, int> > matchs = match.getMatchs();
            ASSERT_EQ(matchs.size(), amnacid.size());
            double score = 0;
            matcher::Score scoring;
            std::vector<aminoacid::Aminoacid> amns;
            int cur_part_id = -1;
            int rev = 0;

            int start = 0;
            if (type == NRP::cycle) {
                for (int i = 0; i < matchs.size(); ++i) {
                    int j = (i - 1 + matchs.size()) % matchs.size();
                    if (matchs[i].first != -1 && (matchs[j].first != matchs[i].first)) {
                        start = i;
                        break;
                    }
                }
            }

            for (int i = 0; i < matchs.size(); ++i) {
                int pos = (start + i) % matchs.size();
                if (matchs[pos].first == -1 &&
                    (i == 0 || matchs[(pos - 1 + matchs.size()) % matchs.size()].first != -1)) {
                    score += scoring.openGap();
                } else if (matchs[pos].first == -1) {
                    score += scoring.continueGap();
                }

                if (cur_part_id != -1 && (matchs[pos].first == -1 || matchs[pos].first != cur_part_id)) {
                    double segscore = 0;
                    if (rev) {
                        std::reverse(amns.begin(), amns.end());
                    }
                    scoring.getScoreForSegment(amns, prediction, cur_part_id, segscore);
                    Segment curseg((pos - amns.size() + matchs.size()) % matchs.size(),
                                   (pos - 1 + matchs.size()) % matchs.size(), 0, 0, segscore);
                    score += scoring.addSegment(curseg);
                    cur_part_id = -1;
                    amns.resize(0);
                    rev = 0;
                }

                if (matchs[pos].first != -1) {
                    nrpsprediction::AAdomainPrediction amn_pred = nrpParts[matchs[pos].first].getAAdomainPrediction()[matchs[pos].second];
                    amns.push_back(amnacid[pos]);
                    cur_part_id = matchs[pos].first;
                    if (matchs[pos].second == 0) {
                        rev = 1;
                    } else {
                        rev = 0;
                    }
                }
            }

            int pos = start;
            if (cur_part_id != -1) {
                double segscore = 0;
                if (rev) {
                    std::reverse(amns.begin(), amns.end());
                }
                scoring.getScoreForSegment(amns, prediction, cur_part_id, segscore);
                Segment curseg((pos - amns.size() + matchs.size()) % matchs.size(),
                               (pos - 1 + matchs.size()) % matchs.size(), 0, 0, segscore);
                score += scoring.addSegment(curseg);
            }

            ASSERT_LE(fabs(score - match.score()), EPS);
        }

        void check_correct_seq(matcher::Matcher::Match match, std::vector< nrpsprediction::OrfPrediction> nrpParts) {
            std::vector<std::vector<int> > post(nrpParts.size());
            for (int i = 0; i < post.size(); ++i) {
                post[i].resize(nrpParts[i].getAAdomainPrediction().size(), -1);
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

        void check_correct_cycle_seq(matcher::Matcher::Match match, std::vector< nrpsprediction::OrfPrediction> nrpParts) {
            std::vector<std::vector<int> > post(nrpParts.size());
            for (int i = 0; i < post.size(); ++i) {
                post[i].resize(nrpParts[i].getAAdomainPrediction().size(), -1);
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


        void check_seq_match_with_nrp(matcher::Matcher::Match match, std::vector< nrpsprediction::OrfPrediction> nrpParts) {
            std::vector<std::pair <int, int> > matchs = match.getMatchs();
            for (int i = 0; i < matchs.size(); ++i) {
                if (matchs[i].first != -1) {
                    nrpsprediction::AAdomainPrediction apr = nrpParts[matchs[i].first].getAAdomainPrediction()[matchs[i].second];
                    ASSERT_TRUE(apr.contain(amnacid[i]));
                }
            }
        }

        nrpsprediction::OrfPrediction getSubPart(int bg, int sz, int delta, double &score) {
            nrpsprediction::OrfPrediction nrps_part("filename", "orf");
            std::vector<aminoacid::Aminoacid> aas;
            int len = amnacid.size();
            int j = 0;
            for (int i = bg; j < sz; i = (i + delta + len)%len, ++j) {
                std::vector<double> prob;
                std::vector<std::string> names;
                for (int g = 0; g < 3; ++g) {
                    prob.push_back(60 + rand()%4 * 10);
                    names.push_back(aminoacid::AminoacidInfo::AMINOACID_NAMES[(rand()%(aminoacid::AminoacidInfo::AMINOACID_CNT - 1))]);
                }

                int right_AA_pos = rand()%3;
                names[right_AA_pos] = amnacid[i].get_name();
                std::sort(prob.rbegin(), prob.rend());
                std::stringstream ss;
                for (int g = 0; g < 3; ++g) {
                    ss << names[g] << "(" << prob[g] << ")";
                    if (g < 2) {
                        ss << ";";
                    }
                }
                using namespace nrpsprediction;
                nrps_part.add_prediction(j + 1, AAdomainPrediction(j + 1,
                                                                   Nrpspredictor2Builder::parse_predictions(ss.str())));
                auto predictions = nrps_part.getAAdomainPrediction();
                aas.push_back(amnacid[i]);
            }

            matcher::Score scoring;
            double curs = 0;
            std::vector<nrpsprediction::OrfPrediction> parts;
            parts.push_back(nrps_part);
            nrpsprediction::BgcPrediction prediction(parts);
            scoring.getScoreForSegment(aas, prediction, 0, curs);
            score += curs;
            
            return nrps_part;
        }
    };


    //check if has segment in NRP then find it.
    TEST_F(ContainNRPsTest, containTrueRandCycleTest) {
        double scr = 0;
        for (int tst = 0; tst < 1000; ++tst) {
            int len = rand()%100 + 2;
            std::shared_ptr<nrp::NRP> nrp = genRandCycleNRP(len);
            int bg = rand()%len, sz = rand()%len + 1;

            int delta = 1;
            if (rand() % 2 == 0) {
                delta = -1;
            }

            int ed = (bg + (sz - 1) * delta + len)%len;

            nrpsprediction::OrfPrediction nrps_part = getSubPart(bg, sz, delta, scr);

            std::vector<nrpsprediction::OrfPrediction> parts;
            parts.push_back(nrps_part);
            matcher::Score score;
            nrpsprediction::BgcPrediction prediction(parts);
            matcher::Matcher matcher1(nrp, &prediction, &score);
            std::vector<Segment> segments = matcher1.matche_seg(0);
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
        matcher::ScoreFullMatch score;
        for (int tst = 0; tst < 1000; ++tst) {
            int len = rand()%20 + 1;
            std::shared_ptr<nrp::NRP> nrp = genRandCycleNRP(len);
            int partlen = rand()%len + 1;

            nrpsprediction::OrfPrediction nrps_part("filename", "orf");
            std::vector<std::vector<aminoacid::Aminoacid>> predict(partlen);

            genRandNrpPart(partlen, nrps_part, predict);

            std::vector<nrpsprediction::OrfPrediction> parts;
            parts.push_back(nrps_part);
            nrpsprediction::BgcPrediction prediction(parts);
            matcher::Matcher matcher1(nrp, &prediction, &score);

            std::vector<Segment> segments = matcher1.matche_seg(0);
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
        matcher::Score score;
        double scr = 0;
        for (int tst = 0; tst < 1000; ++tst) {
            int len = rand()%100 + 2;
            std::shared_ptr<nrp::NRP> nrp = genRandLineNRP(len);

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

            nrpsprediction::OrfPrediction nrps_part = getSubPart(bg, sz, delta, scr);

            std::vector<nrpsprediction::OrfPrediction> parts;
            parts.push_back(nrps_part);
            nrpsprediction::BgcPrediction prediction(parts);
            matcher::Matcher matcher1(nrp, &prediction, &score);

            std::vector<Segment> segments = matcher1.matche_seg(0);
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
        matcher::Score score;
        for (int tst = 0; tst < 1000; ++tst) {
            int len = rand()%20 + 1;
            std::shared_ptr<nrp::NRP> nrp = genRandLineNRP(len);

            int partlen = rand()%len + 1;
            nrpsprediction::OrfPrediction nrps_part ("filename", "orf");
            std::vector<std::vector<aminoacid::Aminoacid>> predict(partlen);

            genRandNrpPart(partlen, nrps_part, predict);

            std::vector<nrpsprediction::OrfPrediction> parts;
            parts.push_back(nrps_part);
            nrpsprediction::BgcPrediction prediction(parts);
            matcher::Matcher matcher1(nrp, &prediction, &score);

            std::vector<Segment> segments = matcher1.matche_seg(0);
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
            std::shared_ptr<nrp::NRP> nrp = genRandLineNRP(len);

            int cntbp = rand()%5 + 1;
            std::vector<int> bps(cntbp);
            for (int i = 0; i < cntbp; ++i) {
                bps[i] = rand()%len;
            }

            bps.push_back(0);
            bps.push_back(len);

            std::sort(bps.begin(), bps.end());

            std::vector<nrpsprediction::OrfPrediction> nrpParts;
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
                nrpsprediction::OrfPrediction nrps_part ("filename", "orf");
                std::vector<std::vector<aminoacid::Aminoacid>> predict(partlen);

                genRandNrpPart(partlen, nrps_part, predict);
                nrpParts.push_back(nrps_part);
            }

            std::random_shuffle(nrpParts.begin(), nrpParts.end());
            nrpsprediction::BgcPrediction nrpsPrediction(nrpParts);

            matcher::ScoreFullMatch score;
            matcher::Matcher matcher(nrp, &nrpsPrediction, &score);
            matcher::Matcher::Match match = matcher.getMatch();

            ASSERT_GE(match.score() - res_score, -EPS);
            check_correct_score(match, nrpParts, NRP::line);
            check_correct_seq(match, nrpParts);
            check_seq_match_with_nrp(match, nrpParts);
        }
    }

    TEST_F(ContainNRPsTest, coverRandCycleTest) {
        matcher::ScoreFullMatch score;
        for (int tst = 0; tst < 1000; ++tst) {
            int len = rand()%20 + 1;
            std::shared_ptr<nrp::NRP> nrp = genRandCycleNRP(len);

            int cntbp = rand()%5 + 1;
            std::vector<int> bps(cntbp);
            for (int i = 0; i < cntbp; ++i) {
                bps[i] = rand()%len;
            }

            std::sort(bps.begin(), bps.end());

            std::vector<nrpsprediction::OrfPrediction> nrpParts;
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
                nrpsprediction::OrfPrediction nrps_part ("filename", "orf");
                std::vector<std::vector<aminoacid::Aminoacid>> predict(partlen);

                genRandNrpPart(partlen, nrps_part, predict);
                nrpParts.push_back(nrps_part);
            }

            std::random_shuffle(nrpParts.begin(), nrpParts.end());
            nrpsprediction::BgcPrediction nrpsPrediction(nrpParts);

            matcher::Matcher matcher(nrp, &nrpsPrediction, &score);
            matcher::Matcher::Match match = matcher.getMatch();

            ASSERT_GE(match.score() - res_score, -EPS);
            check_correct_score(match, nrpParts, NRP::cycle);
            check_correct_cycle_seq(match, nrpParts);
            check_seq_match_with_nrp(match, nrpParts);
        }
    }

    TEST_F(ContainNRPsTest, BranchCycleFakeEnd) {
        typedef aminoacid::Aminoacid Aminoacid;

        //generate branch cycle
        int len = 6;
        std::vector<std::string> strformula(len);
        std::vector<int> first_part = {0, 1};
        std::vector<int> second_part = {2, 3, 4, 5};
        std::vector<Aminoacid> aa1 = {Aminoacid("none"), Aminoacid("val")};
        std::vector<Aminoacid> aa2 = {Aminoacid("trp"), Aminoacid("aad"), Aminoacid("gua"), Aminoacid("orn")};
        std::vector<int> pos1 = first_part, pos2 = first_part;
        std::vector<Aminoacid> aa_sum1 = aa1, aa_sum2 = aa1;
        pos1.insert(pos1.end(), second_part.begin(), second_part.end());
        pos2.insert(pos2.end(), second_part.rbegin(), second_part.rend());
        aa_sum1.insert(aa_sum1.end(), aa2.begin(), aa2.end());
        aa_sum2.insert(aa_sum2.end(), aa2.rbegin(), aa2.rend());


        std::shared_ptr<NRP> nrp1 = std::make_shared<NRPLine>("", strformula, aa_sum1, pos1, "", "");
        std::shared_ptr<NRP> nrp2 = std::make_shared<NRPLine>("", strformula, aa_sum2, pos2, "", "");
        std::shared_ptr<NRP> nrp = std::make_shared<NRPBranch>(nrp1, nrp2, first_part.size());

        nrpsprediction::OrfPrediction nrps_part("", "");

        std::vector<Aminoacid> aas;
        //generate part for prediction
        for (int i = 1; i < aa_sum1.size(); ++i) {
            std::vector<double> prob;
            std::vector<std::string> names;
            for (int g = 0; g < 3; ++g) {
                prob.push_back(60 + rand() % 4 * 10);
                names.push_back(
                        aminoacid::AminoacidInfo::AMINOACID_NAMES[(rand()%(aminoacid::AminoacidInfo::AMINOACID_CNT))]);
            }

            int right_AA_pos = rand() % 3;
            names[right_AA_pos] = aa_sum1[i].get_name();
            std::sort(prob.rbegin(), prob.rend());
            std::stringstream ss;
            for (int g = 0; g < 3; ++g) {
                ss << names[g] << "(" << prob[g] << ")";
                if (g < 2) {
                    ss << ";";
                }
            }
            using namespace nrpsprediction;
            nrps_part.add_prediction(i, AAdomainPrediction(i,
                                                           Nrpspredictor2Builder::parse_predictions(ss.str())));
            aas.push_back(aa_sum1[i]);
        }


        //calculate score
        matcher::Score scoring;
        double score = 0, segscore = 0;

        std::vector<nrpsprediction::OrfPrediction> parts;
        parts.push_back(nrps_part);
        nrpsprediction::BgcPrediction prediction1(parts);
        scoring.getScoreForSegment(aas, prediction1, 0, segscore);
        score += scoring.addSegment(matcher::Segment(1, 5, 0, 0, segscore));

        //match
        nrpsprediction::BgcPrediction prediction(parts);
        matcher::Matcher matcher1(nrp, &prediction, &scoring);
        auto match = matcher1.getMatch();

        //compare score
        ASSERT_LE(fabs(score - match.score()), EPS);
    }
}

int main(int argc, char *argv[]) {
    ::testing::InitGoogleTest(&argc, argv);
    logging::create_logger("");
    return RUN_ALL_TESTS();
}