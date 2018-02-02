#ifndef NRPSMATCHER_NRP_H
#define NRPSMATCHER_NRP_H

#include <vector>
#include "../Aminoacid/Aminoacids.h"
#include "../NRPsPrediction/NRPsPart.h"

namespace nrp {
    class NRP {
    private:
        friend class ContainNRPsTest;
        std::vector <aminoacid::Aminoacids::Aminoacid> aminoacids;
        std::vector <std::string> strformula;
        std::string graph;
        std::string file_name;
    public:
        struct Segment {
            int l;
            int r;
            bool rev;

            Segment() {}
            Segment(int l, int r, bool rev): l(l), r(r), rev(rev) {}
        };

        void parse_fragment_graph(std::string fragment_graph);
        std::vector<Segment> containNRPsPart(nrpsprediction::NRPsPart predict_part);

        int getLen();
        int getInd(int i);
        std::string getFormula(int i);
        aminoacid::Aminoacids::Aminoacid getAminoacid(int i);
        void print();

        std::string getGraphInString();
        std::string get_file_name();
    };
}


#endif //NRPSMATCHER_NRP_H
