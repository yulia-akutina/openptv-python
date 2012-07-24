/*
Implementation the tracking frame buffers. Each frame holds target information
for all cameras, correspondence information and path links information.

Contains functions for dealing with the frame content, such as comparing and
reading targets, correspondences, etc.
*/

#include "tracking_frame_buf.h"

/* Check that target t1 is equal to target t2, i.e. all their fields are equal.
 * 
 * Arguments:
 * target *t1, *t2 - the target structures to compare.
 * 
 * Returns:
 * true if all fields are equal, false otherwise.
*/
int compare_targets(target *t1, target *t2) {
    return (\
        (t1->pnr == t2->pnr) && (t1->x == t2->x) && (t1->y == t2->y) && \
        (t1->n == t2->n) && (t1->nx == t2->nx) && (t1->ny == t2->ny) && \
        (t1->sumg == t2->sumg) && (t1->tnr == t2->tnr));
}

/* Check that two correspondence structs are equal, i.e. all their fields are 
 * equal.
 * 
 * Arguments:
 * corres *c1, *c2 - the target structures to compare.
 * 
 * Returns:
 * true if all fields are equal, false otherwise.
*/
int compare_corres(corres *c1, corres *c2) {
    return ((c1->nr == c2->nr) && (c1->p[0] == c2->p[0]) && \
        (c1->p[1] == c2->p[1]) && (c1->p[2] == c2->p[2]) && \
        (c1->p[3] == c2->p[3]));
}

/* Check that two path info (P) structs are equal, i.e. all their fields are 
 * equal.
 * 
 * Arguments:
 * P *p1, *p2 - the target structures to compare.
 * 
 * Returns:
 * true if all fields are equal, false otherwise.
*/
int compare_path_info(P *p1, P *p2) {
    int equal = 0, iter;
    if (!((p1->prev == p2->prev) && (p1->next == p2->next) && \
        (p1->prio == p2->prio) && (p1->finaldecis == p2->finaldecis) && \
        (p1->inlist == p2->inlist)))
        return 0;

    for (iter = 0; iter < 3; iter++) {
        if (p1->x[iter] != p2->x[iter]) return 0;
    }
    for (iter = 0; iter < POSI; iter++) {
        if (p1->decis[iter] != p2->decis[iter]) return 0;
        if (p1->linkdecis[iter] != p2->linkdecis[iter]) return 0;
    }
    return 1;
}

