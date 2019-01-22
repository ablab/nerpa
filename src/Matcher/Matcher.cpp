//
// Created by olga on 22.01.19.
//

#include "Matcher.h"

nrp::NRP::Match matcher::Matcher::getMatch() const {
    return nrp.isCover(prediction);
}
