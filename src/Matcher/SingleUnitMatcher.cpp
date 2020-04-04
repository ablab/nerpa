//
// Created by olga on 21.07.19.
//

#include "SingleUnitMatcher.h"

namespace matcher {
    matcher::SingleUnitMatcher::SingleUnitMatcher(const std::shared_ptr<nrp::NRP> &nrp,
                                                  const nrpsprediction::BgcPrediction *prediction,
                                                  const matcher::Score *score) : Matcher(nrp, prediction, score) {}

    matcher::SingleUnitMatcher::SingleUnitMatcher() {}

    void matcher::SingleUnitMatcher::matchSingleUnits(matcher::MatcherBase::Match &match, std::vector<bool> &used_pos) const {
        const double INF = 1e9;
        auto short_parts = prediction->getShortOrfs();
        int cnt_big_parts = prediction->getOrfs().size();
        int n = nrp->getLen();
        int m = short_parts.size();
        bool isSwap = false;
        if (n > m) {
            isSwap = true;
            std::swap(n, m);
        }
        std::vector<std::vector<double>> a(n + 1,
                                           std::vector<double>(m + 1, 1));

        for (int i = 1; i <= n; ++i) {
            for (int j = 1; j <= m; ++j) {
                if (!isSwap) {
                    a[i][j] = -score->singleUnitScore(short_parts[j - 1].getAAdomainPrediction()[0],
                                                      nrp->getAminoacid(i - 1));
                    if (used_pos[i - 1] == true) {
                        a[i][j] = 1;
                    }
                } else {
                    a[i][j] = -score->singleUnitScore(short_parts[i - 1].getAAdomainPrediction()[0],
                                                      nrp->getAminoacid(j - 1));
                    if (used_pos[j - 1] == true) {
                        a[i][j] = 1;
                    }
                }
            }
        }

        std::vector<double> u(n + 1), v(m + 1);
        std::vector<int> p(m + 1), way(m + 1);

        for (int i = 1; i <= n; ++i) {
            p[0] = i;
            int j0 = 0;
            std::vector<double> minv(m + 1, INF);
            std::vector<bool> used(m + 1, false);

            do {
                used[j0] = true;
                int i0 = p[j0], delta = INF, j1;
                for (int j=1; j <= m; ++j) {
                    if (!used[j]) {
                        double cur = a[i0][j]-u[i0]-v[j];
                        if (cur < minv[j]) {
                            minv[j] = cur, way[j] = j0;
                        }
                        if (minv[j] < delta) {
                            delta = minv[j], j1 = j;
                        }
                    }
                }

                for (int j = 0; j <= m; ++j) {
                    if (used[j]) {
                        u[p[j]] += delta;
                        v[j] -= delta;
                    } else {
                        minv[j] -= delta;
                    }
                }
                j0 = j1;
            } while (p[j0] != 0);
            do {
                int j1 = way[j0];
                p[j0] = p[j1];
                j0 = j1;
            } while (j0);
        }

        double newScore = match.score();
        for (int j = 1; j <= m; ++j) {
            int npos = p[j];
            int part_pos = j;
            if (npos == 0) continue;
            if (isSwap) {
                std::swap(npos, part_pos);
            }
            if (used_pos[npos - 1] == false) {
                double curs = score->singleUnitScore(short_parts[part_pos - 1].getAAdomainPrediction()[0],
                                                     nrp->getAminoacid(npos - 1));
                if (curs >= 0) {
                    match.match(npos - 1, cnt_big_parts + part_pos - 1, 0);
                    newScore += curs;
                }
            }
        }

        match.setScore(newScore);
    }

    MatcherBase::Match SingleUnitMatcher::updateMatch(const nrpsprediction::BgcPrediction &nrPsPrediction,
                                                      matcher::MatcherBase::Match match, int bg,
                                                      std::vector<Segment> &matched_parts_id) const {
        auto nmatch = setUpdateMatch(nrPsPrediction, match, bg, matched_parts_id);
        std::vector<std::pair<int, int> > part_id_pos = nmatch.getMatchs();
        std::vector<bool> matched_nrp_pos(nrp->getLen(), false);
        for (int i = 0; i < part_id_pos.size(); ++i) {
            if (part_id_pos[i].first > -1) {
                matched_nrp_pos[i] = true;
            }
        }

        matchSingleUnits(nmatch, matched_nrp_pos);
        nmatch.setScore(score->resultScore(nmatch.score(), nrp->getLen(), matched_parts_id, *prediction, *nrp));
        return nmatch;
    }
}