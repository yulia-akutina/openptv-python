/* Utilities for handling sdimple 3d vectors implemented as an array of 3
doubles.
*/

#ifndef VEC_UTILS_H
#define VEC_UTILS_H

#include <math.h>

# define norm(x,y,z) sqrt((x)*(x) + (y)*(y) + (z)*(z))

typedef double pos3d[3];

void init_pos3d(pos3d init);
void copy_pos3d(pos3d dest, pos3d src);
void subst_pos3d(pos3d from, pos3d sub, pos3d output);
double diff_norm_pos3d(pos3d vec1, pos3d vec2);
double dot_pos3d(pos3d vec1, pos3d vec2);

#endif

