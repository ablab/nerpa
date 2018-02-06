#include <fstream>
#include <iostream>
#include <sstream>
#include <algorithm>
#include "NRP.h"

std::string nrp::NRP::getFormula(int i) {
    return strformula[i];
}

std::string nrp::NRP::getGraphInString() {
    return graph;
}

int nrp::NRP::getLen() {
    return aminoacids.size();
}

int nrp::NRP::getInd(int i) {
    return position[i];
}

aminoacid::Aminoacids::Aminoacid nrp::NRP::getAminoacid(int i) {
    return aminoacids[i];
}

void nrp::NRP::print() {
    for (int i = 0; i < (int)aminoacids.size(); ++i) {
        std::cerr << aminoacid::Aminoacids::AMINOACID_NAMES[aminoacids[i]] << " ";
    }
    std::cerr << "\n";
}


std::string nrp::NRP::get_file_name() {
    return file_name;
}

nrp::NRP::Match
nrp::NRP::isCoverLine(std::vector<nrp::NRP::Segment> &segments, nrpsprediction::NRPsPrediction nrPsPrediction,
                      const std::vector<int> &toSmallId, const std::vector<int> &toBigId) {
    int len = this->getLen();
    std::sort(segments.begin(), segments.end());

    std::vector<std::vector<int> > d(len + 1, std::vector<int>((1 << toBigId.size()), -1));
    std::vector<std::vector<std::pair<int, int> > > p(len + 1,
                                                      std::vector<std::pair<int, int> >((1 << toBigId.size()),
                                                                                        std::make_pair(-1, -1)));
    std::vector<std::vector<int> > pa(len + 1, std::vector<int>((1 << toBigId.size()), -1));
    d[0][0] = 0;

    int curseg = 0;

    for (int pos = 0; pos <= len; ++pos) {
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
                    if (((msk >> (toSmallId[segments[curseg].part_id])) & 1) != 0) {
                        ++curseg;
                        continue;
                    }

                    nmsk = msk | (1 << (toSmallId[segments[curseg].part_id]));
                }

                if (d[segments[curseg].r + 1][nmsk] == -1 || d[segments[curseg].r + 1][nmsk] > d[pos][msk]) {
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
        if (mn >= d[len][msk] && d[len][msk] != -1) {
            mn = d[len][msk];
            rmsk = msk;
        }
    }

    Match nrPsMatch(this, nrPsPrediction.getNrpsParts());
    int pos = len;
    while (pos > 0) {
        int nxtp = p[pos][rmsk].first;
        int nxtmsk = p[pos][rmsk].second;
        int seg = pa[pos][rmsk];
        if (seg != -1) {
            int j = segments[seg].r - segments[seg].l;
            if (segments[seg].rev) {
                j = 0;
            }

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
