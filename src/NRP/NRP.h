#ifndef NRPSMATCHER_NRP_H
#define NRPSMATCHER_NRP_H

#include <vector>
#include "../NRPsPrediction/NRPsPrediction.h"
#include "../Aminoacid/Aminoacids.h"
#include "../NRPsPrediction/NRPsPart.h"
#include <iostream>
#include "assert.h"

namespace nrp {
    class NRP;
    class NRPLine;

    class NRP {
    public:
        struct Segment {
            int l;
            int r;
            int part_id;
            bool rev;
            double scor;

            Segment(){}
            Segment(int l, int r, int id, bool rev, double scor): l(l), r(r), part_id(id), rev(rev), scor(scor) {}

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
        std::string extra_info;
    public:

        enum NRPType {cycle, line, branch_cycle};
        NRP() = default;
        NRP(std::string file_name, std::vector<std::string> strformula,
            std::vector <aminoacid::Aminoacids::Aminoacid> aminoacids, std::vector<int> position,
            std::string graph, std::string extra_info):
                file_name(file_name), strformula(strformula), aminoacids(aminoacids), position(position),
                graph(graph), extra_info(extra_info) {
            assert(position.size() == aminoacids.size());
        }

        virtual std::vector<Segment> containNRPsPart(nrpsprediction::NRPsPart predict_part) const = 0;

        virtual int getLen() const;

        virtual int getInd(int i) const;

        virtual std::string getFormula(int i) const;

        virtual aminoacid::Aminoacids::Aminoacid getAminoacid(int i) const;

        virtual void print() const;

        virtual std::string getGraphInString() const;

        virtual std::string get_file_name() const;

        virtual std::string get_extra_info() const;

        virtual NRPType getType() const = 0;

        virtual std::vector<NRPLine> getLines() const = 0;
    };

}


#endif //NRPSMATCHER_NRP_H
