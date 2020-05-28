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
        std::string findAAname() const;
    public:
        enum Configuation{L, D, NA};

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
        std::string get_possible_name() const;
        std::vector<aminoacid::Modification> getModifications() const {
            return modifications;
        }

        int get_id() const;
        bool is_AA() const;
        void addModification(Modification m);

        Configuation getConfiguration() const {
            return configuation;
        }

        void setConfiguration(Configuation conf) {
            configuation = conf;
        }
    private:
        Formula formula;
    public:
        const Formula &getFormula() const;
    private:
        int aa;
        Configuation configuation = NA;
        std::vector<Modification> modifications;
        std::string possible_name = "";
    };
}


#endif //NRPSMATCHER_AMINOACIDS_H
