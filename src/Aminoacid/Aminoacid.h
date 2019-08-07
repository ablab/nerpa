#ifndef NRPSMATCHER_AMINOACIDS_H
#define NRPSMATCHER_AMINOACIDS_H

#include <string>
#include "Formula.h"
#include "Modification.h"
#include "AminoacidInfo.h"
#include "ModificationInfo.h"
#include <vector>

namespace aminoacid {
    class Aminoacid {
    private:
        friend class ModificationTest;
    public:
        explicit Aminoacid(std::string aminoacid_name);
        explicit Aminoacid(Formula formula);
        explicit Aminoacid(int aid);
        Aminoacid();

        bool operator == (const Aminoacid& b) const;
        Formula operator - (const Aminoacid& b) const;
        friend std::ostream& operator << (std::ostream& out, const Aminoacid& a) {
            out << a.get_name();
            for (int i = 0; i < a.modifications.size(); ++i) {
                out << "+" << ModificationInfo::NAMES[a.modifications[i].getId()];
            }
            return out;
        }
        std::string get_name() const;
        void addModification(Modification m);

    private:
        Formula formula;
    public:
        const Formula &getFormula() const;
    private:
        int aa;
        std::vector<Modification> modifications;
    };
}


#endif //NRPSMATCHER_AMINOACIDS_H
