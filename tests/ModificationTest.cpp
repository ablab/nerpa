//
// Created by olga on 03.02.19.
//

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
#include <Matcher/Score/Base/ScoreWithModification.h>
#include <NRPsPrediction/Builders/Nrpspredictor2Builder.h>
#include <memory>

namespace nrp {
    typedef matcher::Segment Segment;
    typedef aminoacid::Modification Modification;
    typedef aminoacid::Aminoacid Aminoacid;
    const double EPS = 1e-4;

    class ModificationTest : public ::testing::Test {
    protected:
        std::vector<aminoacid::Aminoacid> amnacid;
        virtual void SetUp() {
            aminoacid::AminoacidInfo::init("../../resources/aminoacids.tsv", "NRPSPREDICTOR2");
            aminoacid::ModificationInfo::init("../../resources/modifications.tsv");
        }

        void addModificationToAA(Aminoacid &aa) {
            if (aa.get_name() == "asn" || aa.get_name() == "glu") {
                aa.addModification(Modification(aminoacid::ModificationInfo::getIdByNameId("methylation")));
            }
        }

        void genRandNrpPart(int partlen, nrpsprediction::OrfPrediction &nrps_part,
                            std::vector<std::vector<aminoacid::Aminoacid>> &predict) {
            predict.resize(partlen);
            for (int j = 0; j < partlen; ++j) {
                std::vector<double> prob;
                std::vector<std::string> names;
                for (int g = 0; g < 3; ++g) {
                    prob.push_back(rand() % 10 * 10);
                    predict[j].push_back(aminoacid::Aminoacid(rand() % (aminoacid::AminoacidInfo::AMINOACID_CNT - 1)));

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
                amnacid.push_back(aminoacid::Aminoacid(rand() % (aminoacid::AminoacidInfo::AMINOACID_CNT - 1)));
                addModificationToAA(amnacid.back());
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
                amnacid.push_back(aminoacid::Aminoacid(rand() % (aminoacid::AminoacidInfo::AMINOACID_CNT - 1)));

                addModificationToAA(amnacid.back());
            }

            std::random_shuffle(position.begin(), position.end());

            std::shared_ptr<nrp::NRP> res = std::make_shared<nrp::NRPLine>("", strformula, amnacid, position, "", "");
            res->aminoacids = amnacid;
            return res;
        }

        nrpsprediction::OrfPrediction getSubPart(int bg, int sz, int delta, double &score) {
            nrpsprediction::OrfPrediction nrps_part("filename", "orf");
            std::vector<aminoacid::Aminoacid> aas;
            int len = amnacid.size();
            int j = 0;
            for (int i = bg; j < sz; i = (i + delta + len) % len, ++j) {
                std::vector<double> prob;
                std::vector<std::string> names;
                for (int g = 0; g < 3; ++g) {
                    prob.push_back(60 + rand() % 4 * 10);
                    names.push_back(aminoacid::AminoacidInfo::AMINOACID_NAMES[rand() % (aminoacid::AminoacidInfo::AMINOACID_CNT - 1)]);
                }

                int right_AA_pos = rand() % 3;
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

            matcher::Score scoring2;
            matcher::Score* score1 = new matcher::Score;
            matcher::ScoreWithModification scoring(std::unique_ptr<matcher::Score>(std::move(score1)));
            double curs = 0;
            nrpsprediction::BgcPrediction preidction({nrps_part});
            bool found_seg = scoring.getScoreForSegment(aas, preidction, 0, curs);
            double curs2 = 0;
            bool found_seg2 = scoring2.getScoreForSegment(aas, preidction, 0, curs2);
            score += curs;

            return nrps_part;
        }

        void checkSubPart(int bg, int sz, int delta, nrpsprediction::OrfPrediction nrps_part) {
            std::vector<aminoacid::Aminoacid> aas;
            int len = amnacid.size();
            int j = 0;
            for (int i = bg; j < sz; i = (i + delta + len) % len, ++j) {
                aas.push_back(amnacid[i]);
            }

            matcher::Score* score1 = new matcher::Score;
            matcher::ScoreWithModification scoring(std::unique_ptr<matcher::Score>(std::move(score1)));;
            double curs = 0;
            nrpsprediction::BgcPrediction preidction({nrps_part});
            ASSERT_TRUE(scoring.getScoreForSegment(aas, preidction, 0, curs));
        }
    };

    TEST_F(ModificationTest, coverRandLineTest) {
        matcher::Score* score1 = new matcher::Score;
        matcher::ScoreWithModification swm(std::unique_ptr<matcher::Score>(std::move(score1)));
        for (int tst = 0; tst < 1000; ++tst) {
            int len = rand() % 20 + 1;
            std::shared_ptr<NRP> nrp = genRandLineNRP(len);

            int cntbp = rand() % 5 + 1;
            std::vector<int> bps(cntbp);
            for (int i = 0; i < cntbp; ++i) {
                bps[i] = rand() % len;
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
                            checkSubPart(bps[i - 1], (bps[i] - bps[i - 1]), 1, nrpParts.back());
                        } else {
                            nrpParts.push_back(getSubPart(bps[i] - 1, (bps[i] - bps[i - 1]), -1, res_score));
                            checkSubPart(bps[i] - 1, (bps[i] - bps[i - 1]), -1, nrpParts.back());
                        }
                    }
                    res_score -= 1;
                }
            }

            for (int i = 0; i < 10; ++i) {
                int partlen = rand() % 10 + 1;
                nrpsprediction::OrfPrediction nrps_part("filename", "orf");
                std::vector<std::vector<aminoacid::Aminoacid>> predict(partlen);

                genRandNrpPart(partlen, nrps_part, predict);
                nrpParts.push_back(nrps_part);
            }

            std::random_shuffle(nrpParts.begin(), nrpParts.end());
            nrpsprediction::BgcPrediction nrpsPrediction(nrpParts);
            matcher::Matcher matcher(nrp, &nrpsPrediction, &swm);
            matcher::Matcher::Match match = matcher.getMatch();

            ASSERT_GE(match.score() - res_score, -EPS);
        }
    }

    TEST_F(ModificationTest, coverRandCycleTest) {
        matcher::Score* score1 = new matcher::Score;
        matcher::ScoreWithModification swm(std::unique_ptr<matcher::Score>(std::move(score1)));

        for (int tst = 0; tst < 1000; ++tst) {
            int len = rand() % 20 + 1;
            std::shared_ptr<NRP> nrp = genRandCycleNRP(len);

            int cntbp = rand() % 5 + 1;
            std::vector<int> bps(cntbp);
            for (int i = 0; i < cntbp; ++i) {
                bps[i] = rand() % len;
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
                            checkSubPart(bps[pi], (bps[i] - bps[pi] + len) % len, 1, nrpParts.back());
                        } else {
                            nrpParts.push_back(getSubPart((bps[i] - 1 + len) % len, (bps[i] - bps[pi] + len) % len, -1,
                                                          res_score));
                            checkSubPart((bps[i] - 1 + len) % len, (bps[i] - bps[pi] + len) % len, -1, nrpParts.back());
                        }
                    }
                }
            }

            for (int i = 0; i < 10; ++i) {
                int partlen = rand() % 10 + 1;
                nrpsprediction::OrfPrediction nrps_part("filename", "orf");
                std::vector<std::vector<aminoacid::Aminoacid>> predict(partlen);

                genRandNrpPart(partlen, nrps_part, predict);
                nrpParts.push_back(nrps_part);
            }

            std::random_shuffle(nrpParts.begin(), nrpParts.end());
            nrpsprediction::BgcPrediction nrpsPrediction(nrpParts);

            matcher::Matcher matcher(nrp, &nrpsPrediction, &swm);
            matcher::Matcher::Match match = matcher.getMatch();

            ASSERT_GE(match.score() - res_score, -EPS);
        }
    }

    TEST_F(ModificationTest, WithoutSameWithScore) {
        for (int tst = 0; tst < 1000; ++tst) {
            int len = 10;
            amnacid.resize(0);
            for (int i = 0; i < len; ++i) {
                amnacid.push_back(aminoacid::Aminoacid(rand() % (aminoacid::AminoacidInfo::AMINOACID_CNT - 1)));
            }

            nrpsprediction::OrfPrediction nrps_part("filename", "orf");
            int j = 0;
            for (int i = 0; j < len; i += 1, ++j) {
                std::vector<double> prob;
                std::vector<std::string> names;
                for (int g = 0; g < 3; ++g) {
                    prob.push_back(60 + rand() % 4 * 10);
                    names.push_back(aminoacid::AminoacidInfo::AMINOACID_NAMES[rand() % (aminoacid::AminoacidInfo::AMINOACID_CNT - 1)]);
                }

                names[0] = amnacid[i].get_name();
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
            }
            matcher::Score* score1 = new matcher::Score;
            matcher::ScoreWithModification scoring(std::unique_ptr<matcher::Score>(std::move(score1)));

            double curs = 0;
            nrpsprediction::BgcPrediction preidction({nrps_part});
            bool found_seg = scoring.getScoreForSegment(amnacid, preidction, 0, curs);
            matcher::Score scoring2;
            double curs2 = 0;
            bool found_seg2 = scoring2.getScoreForSegment(amnacid, preidction, 0, curs2);
            ASSERT_LE(fabs(curs - curs2), EPS);
        }
    }
}

int main(int argc, char *argv[]) {
    ::testing::InitGoogleTest(&argc, argv);
    logging::create_logger("");
    return RUN_ALL_TESTS();
}
