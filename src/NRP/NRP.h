#ifndef NRPSMATCHER_NRP_H
#define NRPSMATCHER_NRP_H

#include <vector>
#include "../NRPsPrediction/NRPsPrediction.h"
#include "Aminoacid/Aminoacid.h"
#include "../NRPsPrediction/NRPsPart.h"
#include <iostream>
#include "assert.h"

namespace nrp {
    class NRP;
    class NRPLine;

    class NRP {
    protected:
        friend class ContainNRPsTest;
        friend class ModificationTest;
        std::vector <aminoacid::Aminoacid> aminoacids;
        std::vector <int> position;
        std::vector <std::string> strformula;
        std::string graph;
        std::string file_name;
        std::string extra_info;
        size_t len;
    public:

        enum NRPType {cycle, line, branch_cycle};
        NRP() = default;
        NRP(std::string file_name, std::vector<std::string> strformula,
            std::vector <aminoacid::Aminoacid> aminoacids, std::vector<int> position,
            std::string graph, std::string extra_info):
                file_name(file_name), strformula(strformula), aminoacids(aminoacids), position(position),
                graph(graph), extra_info(extra_info) {
            assert(position.size() == aminoacids.size());
            len = aminoacids.size();
        }

        virtual int getLen() const;

        virtual int getFullLen() const;

        virtual int getInd(int i) const;

        virtual std::string getFormula(int i) const;

        virtual aminoacid::Aminoacid getAminoacid(int i) const;

        virtual std::vector<aminoacid::Aminoacid> getAminoacids() const;

        virtual void print() const;

        virtual std::string getGraphInString() const;

        virtual std::string get_file_name() const;

        virtual std::string get_extra_info() const;

        virtual NRPType getType() const = 0;

        virtual std::vector<NRPLine*> getLines() const = 0;

        virtual bool is_valid_seg(int l, int r, int stp) const = 0;

        void deleteAA(int i);
    };

}


#endif //NRPSMATCHER_NRP_H
