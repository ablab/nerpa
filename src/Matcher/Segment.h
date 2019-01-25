//
// Created by olga on 25.01.19.
//

#ifndef NRPSMATCHER_SEGMENT_H
#define NRPSMATCHER_SEGMENT_H

namespace matcher {
    struct Segment {
        int l;
        int r;
        int part_id;
        bool rev;
        double scor;

        Segment() {}

        Segment(int l, int r, int id, bool rev, double scor) : l(l), r(r), part_id(id), rev(rev), scor(scor) {}

        bool operator<(Segment b) {
            return l < b.l;
        }
    };
}


#endif //NRPSMATCHER_SEGMENT_H
