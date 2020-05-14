//
// Created by tag on 24/03/2020.
//

#ifndef NERPA_MONOMERICNRPBUILDER_H
#define NERPA_MONOMERICNRPBUILDER_H

#include "NRPBuilder.h"

namespace nrp {
    class MonomericNRPBuilder : NRPBuilder {
    public:
        static std::shared_ptr<nrp::NRP> build(std::string nrp_id, std::string extra);
    };
}

#endif //NERPA_MONOMERICNRPBUILDER_H
