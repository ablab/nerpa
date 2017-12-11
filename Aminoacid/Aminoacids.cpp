#include <assert.h>
#include "Aminoacids.h"

aminoacid::Aminoacids::Aminoacid aminoacid::Aminoacids::getAminoacid(std::string aminoacid_name) const {
    for (int i = 0; i < AMINOACID_CNT; ++i) {
        if (AMINOACID_NAMES[i] == aminoacid_name) {
            return static_cast<Aminoacid>(i);
        }
    }

    assert(true);
}
