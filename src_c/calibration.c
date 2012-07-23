/* Implementation of calibration methods defined in calibration.h */

#include "calibration.h"

/************************************************
 * compare_exterior() performs a comparison of two Exterior objects.
 *
 * Arguments: 
 * Exterior *e1, *e2 - pointers for the exterior objects to compare.
 *
 * Returns:
 * 1 if equal, 0 if any field differs.
 */
int compare_exterior(Exterior *e1, Exterior *e2) {
    int row, col;
    
    for (row = 0; row < 3; row++)
        for (col = 0; col < 3; col++)
            if (e1->dm[row][col] != e2->dm[row][col])
                return 0;
    
    return ((e1->x0 == e2->x0) && (e1->y0 == e2->y0) && (e1->z0 == e2->z0) \
        && (e1->omega == e2->omega) && (e1->phi == e2->phi) \
        && (e1->kappa == e2->kappa));
}

/************************************************
 * compare_interior() performs a comparison of two Interior objects.
 *
 * Arguments: 
 * Interior *i1, *i2 - pointers for the calibration objects to compare.
 *
 * Returns:
 * 1 if equal, 0 if any field differs.
 */
int compare_interior(Interior *i1, Interior *i2) {
    return ((i1->xh == i2->xh) && (i1->yh == i2->yh) && (i1->cc == i2->cc));
}

/************************************************
 * compare_glass() performs a comparison of two Glass objects.
 *
 * Arguments: 
 * Glass *g1, *g2 - pointers for the Glass objects to compare.
 *
 * Returns:
 * 1 if equal, 0 if any field differs.
 */
int compare_glass(Glass *g1, Glass *g2) {
    return ((g1->vec_x == g2->vec_x) && (g1->vec_y == g2->vec_y) && \
        (g1->vec_z == g2->vec_z));
}

/************************************************
 * compare_addpar() performs a comparison of two ap_52 objects.
 *
 * Arguments: 
 * ap_52 *a1, *a2 - pointers for the calibration objects to compare.
 *
 * Returns:
 * 1 if equal, 0 if any field differs.
 */
int compare_addpar(ap_52 *a1, ap_52 *a2) {
    return ((a1->k1 == a2->k1) && (a1->k2 == a2->k2) && (a1->k3 == a2->k3) && \
        (a1->p1 == a2->p1) && (a1->p2 == a2->p2) && (a1->scx == a2->scx) && \
        (a1->she == a2->she));
}

/************************************************
 * compare_calib() performs a deep comparison of two Calibration objects.
 *
 * Arguments: 
 * Calibration *c1, *c2 - pointers for the calibration objects to compare.
 *
 * Returns:
 * 1 if equal, 0 if any field differs.
 */
int compare_calib(Calibration *c1, Calibration *c2) {
    return (compare_exterior(&(c1->ext_par), &(c2->ext_par)) && \
        compare_interior(&(c1->int_par), &(c2->int_par)) && \
        compare_glass(&(c1->glass_par), &(c2->glass_par)) && \
        compare_addpar(&(c1->added_par), &(c2->added_par)));
}

