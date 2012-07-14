
/*  Unit-test suit for the C core of PyPTV. Uses the Check framework:
    http://check.sourceforge.net/
    
    To run it, type "make check" when in the top C directory, src_c/
*/

#include <check.h>
#include <stdlib.h>

#include "../src_c/ptv.h"

START_TEST(test_allocate_tracking_structs)
{
    // Output buffers:
    target *targs[4][4];
    corres *correspond[4];
    P *path_info[4];
    
    // Dummy things to store in the buffers:
    target t_target;
    corres t_corres;
    P t_path;
    
    int cams = 4;
    int max_targets = 100;
    int frame_ix = 0, cam_ix = 0;
    
    allocate_tracking_structs(targs, correspond, path_info, cams, max_targets);
    fail_unless(trackallocflag == 1);
    
    /* Try to write stuff into the allocated memory and see it doesn't
    segfault.*/
    for (frame_ix = 0; frame_ix < 4; frame_ix++) {
        correspond[frame_ix][42] = t_corres;
        path_info[frame_ix][42] = t_path;
        
        for (cam_ix = 0; cam_ix < cams; cam_ix ++) {
            targs[frame_ix][cam_ix][42] = t_target;
        }
    }
}
END_TEST

Suite* ptv_suite(void) {
    Suite *s = suite_create ("PTV");
    TCase *tc_tw = tcase_create ("Tracking window");
    tcase_add_test(tc_tw, test_allocate_tracking_structs);
    suite_add_tcase (s, tc_tw);
    
    return s;
}

int main(void) {
    int number_failed;
    Suite *s = ptv_suite ();
    SRunner *sr = srunner_create (s);
    srunner_run_all (sr, CK_ENV);
    number_failed = srunner_ntests_failed (sr);
    srunner_free (sr);
    return (number_failed == 0) ? EXIT_SUCCESS : EXIT_FAILURE;
}

