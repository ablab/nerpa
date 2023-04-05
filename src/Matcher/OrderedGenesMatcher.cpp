//
// Created by olga on 13.11.19.
//

#include <Logger/logger.hpp>
#include "OrderedGenesMatcher.h"

template <class T>
std::vector<std::vector<std::vector<T>>> three_dim_vector(int d1, int d2, int d3, T value){
  return std::vector<std::vector<std::vector<T>>> (d1, std::vector<std::vector<T>>(
 						   d2, std::vector<T>(
						   d3, value)));
}

class NRP_iterator_skip_first : public matcher::OrderedGenesMatcher::NRP_iterator {
public:
    explicit NRP_iterator_skip_first(std::shared_ptr<nrp::NRP> nrp) : NRP_iterator(nrp) {}

    int getID(int i) override {
        return i + 1;
    }

    int getLen() override {
        return nrp_->getLen() - 1;
    }
};

class NRP_iterator_skip_last : public matcher::OrderedGenesMatcher::NRP_iterator {
public:
    explicit NRP_iterator_skip_last(std::shared_ptr<nrp::NRP> nrp) : NRP_iterator(nrp) {}

    int getLen() override {
        return nrp_->getLen() - 1;
    }
};

class NRP_iterator_skip_tails : public matcher::OrderedGenesMatcher::NRP_iterator {
public:
    explicit NRP_iterator_skip_tails(std::shared_ptr<nrp::NRP> nrp) : NRP_iterator(nrp) {}

    int getID(int i) override {
        return i + 1;
    }

    int getLen() override {
        return nrp_->getLen() - 2;
    }
};

class NRP_iterator_rev : public matcher::OrderedGenesMatcher::NRP_iterator {
public:
    explicit NRP_iterator_rev(std::shared_ptr<nrp::NRP> nrp):NRP_iterator(nrp){}

    int getID(int i) override {
        return nrp_->getLen() - i - 1;
    }
};


class NRP_iterator_rev_skip_first : public matcher::OrderedGenesMatcher::NRP_iterator {
public:
    explicit NRP_iterator_rev_skip_first(std::shared_ptr<nrp::NRP> nrp) : NRP_iterator(nrp) {}

    int getID(int i) override {
        return nrp_->getLen() - i - 1;
    }

    int getLen() override {
        return nrp_->getLen() - 1;
    }
};

class NRP_iterator_rev_skip_last : public matcher::OrderedGenesMatcher::NRP_iterator {
public:
    explicit NRP_iterator_rev_skip_last(std::shared_ptr<nrp::NRP> nrp) : NRP_iterator(nrp) {}

    int getID(int i) override {
        return nrp_->getLen() - i - 2;
    }

    int getLen() override {
        return nrp_->getLen() - 1;
    }
};

class NRP_iterator_rev_skip_tails : public matcher::OrderedGenesMatcher::NRP_iterator {
public:
    explicit NRP_iterator_rev_skip_tails(std::shared_ptr<nrp::NRP> nrp) : NRP_iterator(nrp) {}

    int getID(int i) override {
        return nrp_->getLen() - i - 2;
    }

    int getLen() override {
        return nrp_->getLen() - 2;
    }
};

class NRP_iterator_cycle : public matcher::OrderedGenesMatcher::NRP_iterator {
private:
    int start_ = 0;
public:
    explicit NRP_iterator_cycle(std::shared_ptr<nrp::NRP> nrp, int start):NRP_iterator(nrp){
        start_ = start;
    }

    int getID(int i) override {
        return (start_ + i) % getLen();
    }
};


class NRP_iterator_cycle_rev : public matcher::OrderedGenesMatcher::NRP_iterator {
private:
    int start_ = 0;
public:
    explicit NRP_iterator_cycle_rev(std::shared_ptr<nrp::NRP> nrp, int start):NRP_iterator(nrp){
        start_ = start;
    }

    int getID(int i) override {
        return (start_ - i + getLen()) % getLen();
    }
};


matcher::MatcherBase::Match matcher::OrderedGenesMatcher::getMatch(const std::shared_ptr<nrp::NRP> nrp,
                                                                   const nrpsprediction::BgcPrediction *prediction,
                                                                   const matcher::Score *score) {

    if (nrp->getType() == nrp::NRP::line) {
        return getLineMatch(true, true, nrp, prediction, score);
    } else if (nrp->getType() == nrp::NRP::cycle) {
        return getCycleMatch(nrp, prediction, score);
    } else {
        return getBranchMatch(nrp, prediction, score);
    }
}


matcher::MatcherBase::Match matcher::OrderedGenesMatcher::getLineMatch(bool can_skip_first, bool can_skip_last,
                                                                       const std::shared_ptr<nrp::NRP> nrp,
                                                                       const nrpsprediction::BgcPrediction *prediction,
                                                                       const matcher::Score *score) const {
    NRP_iterator* nrp_iterator;
    NRP_iterator* nrp_iterator_rev;
    if (can_skip_first && can_skip_last &&
    !nrp->getAminoacid(0).is_AA() && !nrp->getAminoacid(nrp->getLen() - 1).is_AA()) {
        nrp_iterator = new NRP_iterator_skip_tails(nrp);
        nrp_iterator_rev = new NRP_iterator_rev_skip_tails(nrp);
    } else if (can_skip_first && !nrp->getAminoacid(0).is_AA()) {
        nrp_iterator = new NRP_iterator_skip_first(nrp);
        nrp_iterator_rev = new NRP_iterator_rev_skip_first(nrp);
    } else if (can_skip_last && !nrp->getAminoacid(nrp->getLen() - 1).is_AA()) {
        nrp_iterator = new NRP_iterator_skip_last(nrp);
        nrp_iterator_rev = new NRP_iterator_rev_skip_last(nrp);
    } else {
        nrp_iterator = new NRP_iterator(nrp);
        nrp_iterator_rev = new NRP_iterator_rev(nrp);
    }

    matcher::MatcherBase::Match matche1 = getSimpleMatch(nrp_iterator, prediction, score);
//    matcher::MatcherBase::Match matche2 = getSimpleMatch(nrp_iterator_rev, prediction, score);

    delete nrp_iterator;
    delete nrp_iterator_rev;

    return matche1;
//    if (matche1.score() > matche2.score()) {
//        return matche1;
//    }
//    return matche2;
}

matcher::MatcherBase::Match matcher::OrderedGenesMatcher::getCycleMatch(const std::shared_ptr<nrp::NRP> nrp,
                                                                        const nrpsprediction::BgcPrediction *prediction,
                                                                        const matcher::Score *score) const {
    matcher::MatcherBase::Match resMatche(nrp, prediction->getOrfs(), score->minScore(nrp->getLen()), score);

    for (int i = 0; i < nrp->getLen(); ++i) {
        NRP_iterator* nrp_iterator = new NRP_iterator_cycle(nrp, i);
        NRP_iterator* nrp_iterator_rev = new NRP_iterator_cycle_rev(nrp, i);

        matcher::MatcherBase::Match matche1 = getSimpleMatch(nrp_iterator, prediction, score);
//        matcher::MatcherBase::Match matche2 = getSimpleMatch(nrp_iterator_rev, prediction, score);

        delete nrp_iterator;
        delete nrp_iterator_rev;

        if (matche1.score() > resMatche.score()) {
            resMatche = matche1;
        }

//        if (matche2.score() > resMatche.score()) {
//            resMatche = matche2;
//        }
    }

    return resMatche;
}

matcher::MatcherBase::Match matcher::OrderedGenesMatcher::getBranchMatch(const std::shared_ptr<nrp::NRP> nrp,
                                                                         const nrpsprediction::BgcPrediction *prediction,
                                                                         const matcher::Score *score) const {

    auto matche1 = getLineMatch(true, false, nrp->getLines()[0], prediction, score);
    auto matche2 = getLineMatch(true, false, nrp->getLines()[1], prediction, score);

    if (matche1.score() > matche2.score()) {
        return matche1;
    }
    return matche2;
}

matcher::MatcherBase::Match matcher::OrderedGenesMatcher::getSimpleMatch(matcher::OrderedGenesMatcher::NRP_iterator *nrp_iterator,
                                                                         const nrpsprediction::BgcPrediction *prediction,
                                                                         const matcher::Score *score) const {
    const auto &nrp_parts = prediction->getOrfs();
    std::vector<int> part_id;
    std::vector<int> pos_id;
    size_t plen = 0;

    for (int pid = 0; pid < nrp_parts.size(); ++pid) {
        size_t cur_len = nrp_parts[pid].getAAdomainPrediction().size();
        for (int i = 0; i < cur_len; ++i) {
            part_id.push_back(pid);
            pos_id.push_back(i);
        }

        plen += cur_len;
    }

    auto get_aapred = [&nrp_parts, &part_id, &pos_id](int ppos) {
        return nrp_parts[part_id[ppos]].getAAdomainPrediction()[pos_id[ppos]];
    };

    auto is_rep_aa = [&nrp_parts, &part_id, &pos_id](int ppos) {
        return nrp_parts[part_id[ppos]].getAAdomainPrediction()[pos_id[ppos]].is_rep();
    };

    auto is_rep_orf = [&nrp_parts, &part_id, &pos_id](int ppos) {
        return nrp_parts[part_id[ppos]].is_rep();
    };

    auto orf_len = [&nrp_parts, &part_id, &pos_id](int ppos) {
        return nrp_parts[part_id[ppos]].getAAdomainPrediction().size();
    };

    auto nrplen = size_t(nrp_iterator->getLen());

    enum match_type {Match, Insertion, Deletion, NA};

    const int max_number_aa_reps = 3;
    const int max_number_orf_reps = 3;

										   
    //dp[x][y][r] is score for aligning nrp[:x] with pred[:y] with the condition that the last orf was repeated r times
    auto dp = three_dim_vector<double>(nrplen + 1, plen + 1, max_number_orf_reps + 1, score->minScore(nrplen));
    //pp[x][y][r] is the ancestor of state (x, y, r) in the optimal route
    auto pp = three_dim_vector<std::tuple<int, int, int, match_type>> (nrplen + 1, plen + 1, max_number_orf_reps + 1,
								      std::make_tuple(0, 0, 0, NA));
    auto recalc = [&](int curx, int cury, int curr,
		      int ancx, int ancy, int ancr,
		      double edge_weight, match_type mt){
		    if (dp[curx][cury][curr] < dp[ancx][ancy][ancr] + edge_weight){
		      dp[curx][cury][curr] = dp[ancx][ancy][ancr] + edge_weight;
		      pp[curx][cury][curr] = std::make_tuple(ancx, ancy, ancr, mt);
		    }
		  };
    
    dp[0][0][0] = 0;
    //in the initial version it was not checked that the new value of dp is greater than the old one. I hope, that is always the case.
    for (int npos = 1; npos <= nrplen; ++npos) 
      recalc(npos,     0, 0,
	     npos - 1, 0, 0,
	     score->InsertionScore(), Insertion);

    //in the initial version it was not checked that the new value of dp is greater than the old one. I hope, that is always the case.
    for (int ppos = 1; ppos <= plen; ++ppos) 
      recalc(0, ppos,     0,
	     0, ppos - 1, 0,
	     score->DeletionScore(), Deletion);

    for (int npos = 1; npos <= nrplen; ++npos) {
        for (int ppos = 1; ppos <= plen; ++ppos) {
	  double match_score =  score->aaScore(get_aapred(ppos - 1), nrp_iterator->getAA(npos - 1));
	  bool new_orf = (ppos == 1) || (part_id[ppos - 1] != part_id[ppos - 2]); //current position in predicted nrp is the beginning of new orf

	  for (int orf_reps = 0; orf_reps < max_number_orf_reps; ++orf_reps)
	    recalc(npos, ppos, new_orf ? 0 : orf_reps,
		   npos - 1, ppos - 1, orf_reps,
		   match_score, Match);

            if (is_rep_aa(ppos - 1)){
	      int max_reps = std::min(max_number_aa_reps, npos);
	      std::vector <double> rep_score(max_reps + 1);
	      rep_score[0] = 0;
	      for (int nreps = 1; nreps <= max_reps; ++nreps)
		rep_score[nreps] = rep_score[nreps - 1] + score->aaScore(get_aapred(ppos - 1),
									 nrp_iterator->getAA(npos - nreps));

	      for(int nreps = 2; nreps <= max_reps; ++nreps)
		for (int orf_reps = 0; orf_reps < max_number_orf_reps; ++orf_reps)
		  recalc(npos, ppos, new_orf ? 0 : orf_reps,
			 npos - nreps, ppos - 1, orf_reps,
			 rep_score[nreps], Match);
	    }

	    if (is_rep_orf(ppos - 1) && pos_id[ppos - 1] == 0 && orf_len(ppos - 1) > 1) //if orf_len == 1, pumping a single domain will do the same job
	      for (int orf_reps = 1; orf_reps < max_number_orf_reps; ++orf_reps)
		recalc(npos, ppos, orf_reps,
		       npos - 1, ppos + orf_len(ppos - 1) - 1, orf_reps - 1,
		       match_score, Match);

	    for (int orf_reps = 0; orf_reps < max_number_orf_reps; ++orf_reps)
	      recalc(npos, ppos, new_orf ? 0 : orf_reps,
		     npos, ppos - 1, orf_reps,
		     score->DeletionScore(), Deletion);

	    for (int orf_reps = 0; orf_reps < max_number_orf_reps; ++orf_reps)
              recalc(npos, ppos, orf_reps,
		     npos - 1, ppos, orf_reps,
                     score->InsertionScore(), Insertion);

	}
    }
    
    std::vector<double>& dp_last = dp[nrplen][plen];
    int reps_last = std::distance(dp_last.begin(),
				   max_element(dp_last.begin(), dp_last.end()));

    int curx = int(nrplen), cury = int(plen), curr = reps_last; 
    //backtracking the optimal match
    matcher::MatcherBase::Match nrPsMatch(nrp_iterator->getNRP(), nrp_parts, dp[curx][cury][curr], score);

//    std::cout << nrp_iterator->getNRP()->get_file_name() << "\n";

    while (curx != 0 || cury != 0) {
      auto anc = pp[curx][cury][curr];
      int ancx = std::get<0> (anc), ancy = std::get<1> (anc), ancr = std::get<2> (anc);
      match_type cur_match = std::get<3>(anc);
        if (cur_match == Match) {
	  for (int xpos = curx - 1; xpos >= ancx; --xpos){
	    nrPsMatch.match_align(nrp_iterator->getID(xpos), part_id[cury - 1], pos_id[cury - 1]);
	    nrPsMatch.match(nrp_iterator->getID(xpos), part_id[cury - 1], pos_id[cury - 1]);
	  }
//            std::cout << "M " << curx << " " << nrp_iterator->getAA(curx - 1).get_name() << "  -- " << cury << " "
//            << get_aapred(cury - 1).getAAPrediction()[0].aminoacid.get_name() << "\n";
        } else if (cur_match == Deletion) {
            nrPsMatch.match_align(-1, part_id[cury - 1], pos_id[cury - 1]);
//            std::cout << "D " << curx << " none -- " << cury << " "
//                                                    << get_aapred(cury - 1).getAAPrediction()[0].aminoacid.get_name() << "\n";
        } else if (cur_match == Insertion) {
            nrPsMatch.match_align(nrp_iterator->getID(curx - 1), -1, -1);
//            std::cout << "I " << curx << " " << nrp_iterator->getAA(curx - 1).get_name() << "  -- " << cury << " none \n";
        }

	curx = ancx, cury = ancy, curr = ancr;
    }

//    std::cout << "\n";
    nrPsMatch.setScore(score->resultScore(dp[nrplen][plen][reps_last], nrplen));
    return nrPsMatch;
    
}
