//
// Created by olga on 22.01.19.
//

#include <Logger/logger.hpp>
#include "Score.h"
#include <Aminoacid/MonomerInfo.h>
#include <math.h>
#include <fstream>

bool matcher::Score::getScoreForSegment(const std::vector<aminoacid::Aminoacid>& amns,
                                        const nrpsprediction::BgcPrediction& prediction, int part_id,
                                        double& score) const {
    if (baseScore != nullptr) {
        return baseScore->getScoreForSegment(amns, prediction, part_id, score);
    } else {
        nrpsprediction::OrfPrediction part = prediction.getOrfs()[part_id];
        std::vector<nrpsprediction::AAdomainPrediction> aminoacid_predictions = part.getAAdomainPrediction();
        int cnt_mismatch = 0;
        int g = 0;
        double segscor = 0;
        for (int j = 0; j < (int) aminoacid_predictions.size() && cnt_mismatch < 2; ++j) {
            if (!aminoacid_predictions[j].contain(amns[j])) {
                cnt_mismatch += 1;
            }
            double cur_sc = aaScore(aminoacid_predictions[j], amns[j]);

            segscor += cur_sc;
        }

        if (cnt_mismatch == 0 || (cnt_mismatch == 1 && aminoacid_predictions.size() > 4)) {
            score = segscor;
            return true;
        }

        return false;
    }
}

double matcher::Score::aaScore(const nrpsprediction::AAdomainPrediction &apred,
                               const aminoacid::Aminoacid &aminoacid) const {
    if (baseScore != nullptr) {
        return baseScore->aaScore(apred, aminoacid);
    } else {
        std::pair<int, int> position = apred.getAmnAcidPos(aminoacid);
        nrpsprediction::AAdomainPrediction::AminoacidProb prob = apred.getAminoacid(aminoacid);
        double score = 0;
        if (getScore(aminoacid, aminoacid, prob, position, score)) {
            return score;
        } else {
            return Mismatch(aminoacid, apred);
        }
    }
}

matcher::Score::Score(std::unique_ptr<Score> base) : Score() {
    baseScore = std::move(base);
}

std::pair<double, aminoacid::Aminoacid>
matcher::Score::getTheBestAAInPred(const nrpsprediction::AAdomainPrediction &apred,
                                   const aminoacid::Aminoacid &aminoacid,
                                   nrpsprediction::AAdomainPrediction::AminoacidProb &probRes,
                                   std::pair<int, int> &posRes) const {
    if (baseScore != nullptr) {
        return baseScore->getTheBestAAInPred(apred, aminoacid, probRes, posRes);
    } else {
        double score = aaScore(apred, aminoacid);
        auto curpos = apred.getAmnAcidPos(aminoacid);
        std::swap(posRes, curpos);
        probRes = apred.getAminoacid(aminoacid);
        return std::pair<double, aminoacid::Aminoacid>(score, probRes.aminoacid);
    }
}

double getModificationScore(const aminoacid::Aminoacid &nrpAA,
        const std::vector<aminoacid::Modification> &mods) {

    std::vector<int> pred_mods_flags(aminoacid::ModificationInfo::MODIFICATION_CNT, 0);
    for (auto &m : mods) {
        pred_mods_flags[m.getId()] = 1;
    }

    auto mds = nrpAA.getModifications();
    std::vector<int> nrp_mods_flags(aminoacid::ModificationInfo::MODIFICATION_CNT, 0);
    for (auto &m : mds) {
        nrp_mods_flags[m.getId()] = 1;
    }

    double res = 0;
    size_t DL_id = aminoacid::ModificationInfo::getIdByNameId("@D");
    for (int i = 0; i != aminoacid::ModificationInfo::MODIFICATION_CNT; ++i) {
        if (i != DL_id) {
            res += aminoacid::ModificationInfo::getCoefficientById(i, pred_mods_flags[i], nrp_mods_flags[i]);
        }
    }
    return res;
//    bool pred_has_met = false;
//    bool nrp_has_met = false;
//    for (int i = 0; i < mods.size(); ++i) {
//        if (mods[i].getId() == aminoacid::ModificationInfo::getIdByNameId("methylation")) {
//            pred_has_met = true;
//        }
//    }
//
//    auto mds = nrpAA.getModifications();
//    for (int i = 0; i < mds.size(); ++i) {
//        if (mds[i].getId() == aminoacid::ModificationInfo::getIdByNameId("methylation")) {
//            nrp_has_met = true;
//        }
//    }
//    if (!nrp_has_met && !pred_has_met) {
//        return 0.1;
//    }
//    if (nrp_has_met && pred_has_met) {
//        return 2.2;
//    }
//    if (pred_has_met && !nrp_has_met) {
//        return -2.2;
//    }
//    if (!pred_has_met && nrp_has_met) {
//        return -2.2;
//    }
}

double matcher::Score::probGenAA(const aminoacid::Aminoacid &nrpAA) const {
    return aminoacid::MonomerInfo::getLogP(nrpAA.get_id());

    std::vector<double> l_list;
    for (int i = 0; i < ProbGenCorrect.size(); ++i) {
        double cur_score = ProbGetScore[i] + ProbGenCorrect[i] + aminoacid::AminoacidInfo::LogP[nrpAA.get_id()];
        //std::cout <<"AA" << nrpAA.get_name() << " " << nrpAA.get_id() << " " << aminoacid::AminoacidInfo::LogP[nrpAA.get_id()] << "\n";
        l_list.push_back(cur_score);
        cur_score = ProbGetScore[i]  + ProbGenIncorrect[i] + aminoacid::MonomerInfo::getLogP(nrpAA.get_id());
        //std::cout << "Mono" << nrpAA.get_name() << " " << nrpAA.get_id() << " " <<  aminoacid::MonomerInfo::getLogP(nrpAA.get_id()) << "\n";

        l_list.push_back(cur_score);
    }

    double exp_val = 0;
    for (int i = 0; i < l_list.size(); ++i) {
        exp_val += exp(l_list[i]);
    }

    //std::cout << exp_val << "\n";
    //std::cout << log(exp_val) << "\n";

    return log(exp_val);
}


double getLDScore(const aminoacid::Aminoacid &nrpAA, const aminoacid::Aminoacid &predAA) {
    auto nrpConf = nrpAA.getConfiguration();
    auto predConf = predAA.getConfiguration();

    if (nrpConf == aminoacid::Aminoacid::NA || predConf == aminoacid::Aminoacid::NA) {
        return 0;
    }

    size_t DL_id = aminoacid::ModificationInfo::getIdByNameId("@D");
    return aminoacid::ModificationInfo::getCoefficientById(DL_id, predConf == aminoacid::Aminoacid::D, nrpConf == aminoacid::Aminoacid::D);
}

bool matcher::Score::getScore(const aminoacid::Aminoacid &nrpAA, const aminoacid::Aminoacid &predAA,
                                const nrpsprediction::AAdomainPrediction::AminoacidProb &prob,
                                const std::pair<int, int> &pos,
                                double& score) const {
    if (baseScore != nullptr) {
        return baseScore->getScore(nrpAA, predAA, prob, pos, score);
    } else {
        if (pos.first != 0 or prob.prob < 68.) {
            return false;
        } else {
            double md_score = getModificationScore(nrpAA, prob.modificatins);
            double dl_score = getLDScore(nrpAA, prob.aminoacid);
            int ind = std::min(int(10 - prob.prob/10), 4);
            //std::cout << nrpAA.get_name() << " " << ind << " " << Score::ProbGenCorrect[ind] << "\n";
            score = Score::ProbGenCorrect[ind];
            score += md_score;
            score += dl_score;
            score -= probGenAA(nrpAA);
            //std::cout << probGenAA(nrpAA) << "\n";
        }
    }

    return true;
}

double matcher::Score::Mismatch(const aminoacid::Aminoacid &structure_aa,
                                const nrpsprediction::AAdomainPrediction &aa_prediction) const {
    if (baseScore != nullptr) {
        return baseScore->Mismatch(structure_aa, aa_prediction);
    } else {
        double md_score = getModificationScore(structure_aa, aa_prediction.get_modifications());
        double dl_score = getLDScore(structure_aa, aa_prediction.getAAPrediction()[0].aminoacid);


        double cur_prob = 60;
        if (!aa_prediction.getAAPrediction().empty()) {
            cur_prob = aa_prediction.getAAPrediction()[0].prob;
        }
        int ind = std::min(int(10 - cur_prob/10), 4);
        double score = Score::ProbGenIncorrect[ind] + aminoacid::MonomerInfo::getLogP(structure_aa.get_id());
        score -= probGenAA(structure_aa);
        return (score + md_score + dl_score);
    }
}

matcher::Score::Score(const std::string &config) {
    std::ifstream in(config);
    std::string tmp;
    in >> tmp >> insertion >> tmp >> deletion;

    auto read_line_to_vec = [&in] (size_t n) -> std::vector<double> {
        std::vector<double> res;
        double w;
        std::string tmp;
        in >> tmp;
        for (size_t i = 0; i != n; ++i) {
            in >> w;
            res.push_back(w);
        }
        return res;
    };

    size_t n_params = 5;
    ProbGenCorrect = read_line_to_vec(n_params);
    ProbGenIncorrect = read_line_to_vec(n_params);
    ProbGetScore = read_line_to_vec(n_params);
}
