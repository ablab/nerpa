#ifndef NRPSMATCHER_NRP_H
#define NRPSMATCHER_NRP_H

#include <vector>
#include "../NRPsPrediction/NRPsPrediction.h"
#include "../Aminoacid/Aminoacids.h"
#include "../NRPsPrediction/NRPsPart.h"
#include <iostream>

namespace nrp {
    class NRP;

    class NRP {
    public:
        class Match {
        private:
            nrp::NRP* nrp;
            std::vector<nrpsprediction::NRPsPart> nrpParts;
            std::vector<int> parts_id;
            std::vector<int> parts_pos;
        public:
            Match() = default;

            Match(nrp::NRP* nrp, std::vector<nrpsprediction::NRPsPart> nrpParts): nrp(nrp), nrpParts(nrpParts) {
                parts_id.resize(nrp->getLen(), -1);
                parts_pos.resize(nrp->getLen(), -1);
            }

            void match(int pos, int part_id, int part_pos);
            void print(std::ofstream& out);
            int score();
            std::vector<std::pair<int, int> > getMatchs();

            bool operator < (Match b);
        };

        struct Segment {
            int l;
            int r;
            int part_id;
            bool rev;

            Segment(){}
            Segment(int l, int r, int id, bool rev): l(l), r(r), part_id(id), rev(rev) {}

            bool operator < (Segment b) {
                return l < b.l;
            }
        };

    protected:
        friend class ContainNRPsTest;
        std::vector <aminoacid::Aminoacids::Aminoacid> aminoacids;
        std::vector <int> position;
        std::vector <std::string> strformula;
        std::string graph;
        std::string file_name;

        virtual Match isCoverLine(std::vector<Segment>& segments, nrpsprediction::NRPsPrediction nrPsPrediction,
                                         const std::vector<int>& toSmallId, const std::vector<int>& toBigId);
    public:
        NRP() = default;
        NRP(std::string file_name, std::vector<std::string> strformula,
            std::vector <aminoacid::Aminoacids::Aminoacid> aminoacids, std::vector<int> position,
            std::string graph): file_name(file_name), strformula(strformula), aminoacids(aminoacids), position(position), graph(graph) {}

        virtual std::vector<Segment> containNRPsPart(nrpsprediction::NRPsPart predict_part) = 0;

        int getLen();
        int getInd(int i);
        std::string getFormula(int i);
        aminoacid::Aminoacids::Aminoacid getAminoacid(int i);
        void print();

        std::string getGraphInString();
        std::string get_file_name();

        virtual Match isCover(nrpsprediction::NRPsPrediction nrPsPrediction) = 0;
    };

}


#endif //NRPSMATCHER_NRP_H
