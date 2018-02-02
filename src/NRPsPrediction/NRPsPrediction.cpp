#include <fstream>
#include <sstream>
#include <iostream>
#include <algorithm>
#include <NRPsMatch/NRPsMatch.h>
#include "NRPsPrediction.h"

namespace nrpsprediction {
    void nrpsprediction::NRPsPrediction::read_file(std::string file_name) {
        std::ifstream in(file_name);
        std::string s;
        while (getline(in, s)) {
            std::stringstream ss(s);
            std::string orf_name;
            std::string predict_aminoacid;
            std::string predict_aminoacids;

            ss >> orf_name >> predict_aminoacid >> predict_aminoacids;
            std::pair <std::string, int> orf_name_num = get_orf_name_and_order(orf_name);

            if (!nrpparts.empty() && nrpparts[nrpparts.size() - 1].get_orf_name() == orf_name_num.first) {
                nrpparts[nrpparts.size() - 1].add_prediction(orf_name_num.second, predict_aminoacids);
            } else {
                if (nrpparts.size() > 0 && nrpparts[nrpparts.size() - 1].getAminoacidsPrediction().size() < 3) {
                    nrpparts.pop_back();
                }
                nrpparts.push_back(NRPsPart(file_name, orf_name_num.first, orf_name_num.second, predict_aminoacids));
            }
        }
        if (nrpparts.size() > 0 && nrpparts[nrpparts.size() - 1].getAminoacidsPrediction().size() < 3) {
            nrpparts.pop_back();
        }


        in.close();
    }

    std::pair<std::string, int> NRPsPrediction::get_orf_name_and_order(std::string orf) {
        std::stringstream ss(orf);
        std::string prefix;
        std::string orfname;
        std::string tail;
        getline(ss, prefix, '_');
        getline(ss, orfname, '_');
        getline(ss, tail, '_');

        std::stringstream ss2(tail.substr(1));
        int pos;
        ss2 >> pos;

        return std::make_pair(orfname, pos);
    }

    std::vector<NRPsPart> NRPsPrediction::getNrpsParts() {
        return nrpparts;
    }

    nrpsmatch::NRPsMatch NRPsPrediction::isCover(nrp::NRP nrp) {
        std::vector<Segment> segments;
        std::vector<int> toSmallId(nrpparts.size(), -1);
        std::vector<int> toBigId;
        int len = nrp.getLen();

        for (int i = 0; i < nrpparts.size(); ++i) {
            std::vector<nrp::NRP::Segment> part_seg = nrp.containNRPsPart(nrpparts[i]);
            for (int j = 0; j < part_seg.size(); ++j) {
                segments.push_back(Segment(part_seg[j].l, part_seg[j].r, i, part_seg[j].rev));
            }

            if (part_seg.size() >= 2) {
                toSmallId[i] = toBigId.size();
                toBigId.push_back(i);
            }
        }

        for (int bg = 0; bg < len; ++bg) {
            std::vector<Segment> tmpSeg;
            for (int i = 0; i < segments.size(); ++i) {
                tmpSeg.push_back(Segment(segments[i].l - bg, segments[i].r - bg, segments[i].part_id, segments[i].rev));
                if (tmpSeg[tmpSeg.size() - 1].l < 0) {
                    tmpSeg[tmpSeg.size() - 1].l += len;
                    tmpSeg[tmpSeg.size() - 1].r += len;
                }
            }

            nrpsmatch::NRPsMatch curMatch = isCoverLine(tmpSeg, nrp, toSmallId, toBigId);

            if (curMatch.score() > len - 2) {
                std::vector<std::pair<int, int> > matchs = curMatch.getMatchs();
                nrpsmatch::NRPsMatch resMatchs(nrp, this->nrpparts);
                for (int i = 0; i < matchs.size(); ++i) {
                    if (matchs[i].first != -1) {
                        resMatchs.match((i + bg + len) % len, matchs[i].first, matchs[i].second);
                    }
                }

                return resMatchs;
            }
        }


        nrpsmatch::NRPsMatch resMatch(nrp, this->nrpparts);
        return resMatch;
    }

    nrpsmatch::NRPsMatch NRPsPrediction::isCoverLine(std::vector<NRPsPrediction::Segment> &segments, nrp::NRP nrp,
                                     const std::vector<int> &toSmallId, const std::vector<int> &toBigId) {
        int len = nrp.getLen();
        std::sort(segments.begin(), segments.end());

        std::vector<std::vector<int> > d(len + 1, std::vector<int>((1 << toBigId.size()), -1));
        std::vector<std::vector<std::pair<int, int> > > p(len + 1, std::vector<std::pair<int, int> >((1 << toBigId.size()), std::make_pair(-1, -1)));
        std::vector<std::vector<int > > pa(len + 1, std::vector<int>((1 << toBigId.size()), -1));
        d[0][0] = 0;

        int curseg = 0;

        for (int pos = 0; pos < len; ++pos) {
            int lstseg = curseg;
            for (int msk = 0; msk < (1 << (toBigId.size())); ++msk) {
                if (pos != 0 && (d[pos][msk] == -1 || d[pos][msk] > d[pos - 1][msk] + 1)) {
                    if (d[pos - 1][msk] != -1) {
                        d[pos][msk] = d[pos - 1][msk] + 1;
                        p[pos][msk].first = pos - 1;
                        p[pos][msk].second = msk;
                        pa[pos][msk] = -1;
                    }
                }
                if (d[pos][msk] == -1) {
                    continue;
                }
                curseg = lstseg;

                while (curseg < segments.size() && segments[curseg].l == pos) {
                    if (segments[curseg].r >= len) {
                        ++curseg;
                        continue;
                    }
                    int nmsk = msk;

                    if (toSmallId[segments[curseg].part_id] != -1) {
                        if (((msk >> (toSmallId[segments[curseg].part_id]))&1) != 0) {
                            ++curseg;
                            continue;
                        }

                        nmsk = msk | (1 << (toSmallId[segments[curseg].part_id]));
                    }

                    if (d[segments[curseg].r + 1][nmsk] == -1) {
                        d[segments[curseg].r + 1][nmsk] = d[pos][msk];
                        p[segments[curseg].r + 1][nmsk].first = pos;
                        p[segments[curseg].r + 1][nmsk].second = msk;
                        pa[segments[curseg].r + 1][nmsk] = curseg;
                    } else if (d[segments[curseg].r + 1][nmsk] > d[pos][msk]) {
                        d[segments[curseg].r + 1][nmsk] = d[pos][msk];
                        p[segments[curseg].r + 1][nmsk].first = pos;
                        p[segments[curseg].r + 1][nmsk].second = msk;
                        pa[segments[curseg].r + 1][nmsk] = curseg;
                    }

                    ++curseg;
                }

            }
        }

        int mn = len;
        int rmsk = 0;
        for (int msk = 0; msk < (1 << (toBigId.size())); ++msk) {
            mn = std::min(mn, d[len][msk]);
            rmsk = msk;
        }

        nrpsmatch::NRPsMatch nrPsMatch(nrp, this->nrpparts);
        int pos = len;
        while (pos > 0) {
            int nxtp = p[pos][rmsk].first;
            int nxtmsk = p[pos][rmsk].second;
            int seg = pa[pos][rmsk];
            int j = segments[seg].r - segments[seg].l;
            if (segments[seg].rev) {
                j = 0;
            }
            if (seg != -1) {
                int curp = pos - 1;
                while (curp >= nxtp) {
                    nrPsMatch.match(curp, segments[seg].part_id, j);
                    if (segments[seg].rev) {
                        ++j;
                    } else {
                        --j;
                    }
                    --curp;
                }
            }

            pos = nxtp;
            rmsk = nxtmsk;
        }

        return nrPsMatch;
    }
}
