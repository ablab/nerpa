#ifndef NRPSMATCHER_MODIFICATIONINFO_H
#define NRPSMATCHER_MODIFICATIONINFO_H

namespace aminoacid {
    class ModificationInfo {
    public:
        static int MODIFICATION_CNT;
        static std::vector<std::string> NAMES;
        static std::vector<Formula> FORMULS;

        static void init(std::string filename);
    };
}


#endif //NRPSMATCHER_MODIFICATIONINFO_H
