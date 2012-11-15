/* implemnentation of vec_utils.h */

#include "vec_utils.h"

void init_pos3d(pos3d init) {
    int ix;
    for (ix = 0; ix < 3; ix ++) init[ix] = -999;
}

void copy_pos3d(pos3d dest, pos3d src) {
    int ix;
    for (ix = 0; ix < 3; ix ++) dest[ix] = src[ix];
}

void subst_pos3d(pos3d from, pos3d sub, pos3d output) {
    int ix;
    for (ix = 0; ix < 3; ix ++) output[ix] = from[ix] - sub[ix];
}

double diff_norm_pos3d(pos3d vec1, pos3d vec2) {
    return norm(vec1[0] - vec2[0], vec1[1] - vec2[1], vec1[2] - vec2[2]);
}

double dot_pos3d(pos3d vec1, pos3d vec2) {
    int ix;
    double sum = 0;
    
    for (ix = 0; ix < 3; ix ++) sum += vec1[ix]*vec2[ix];
    return sum;
}

