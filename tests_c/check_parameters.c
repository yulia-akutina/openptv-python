/* Unit tests for reading and writing parameter files. */

#include <check.h>
#include <stdlib.h>
#include <stdio.h>
#include "../src_c/parameters.h"

START_TEST(test_read_sequence_par)
{
    int cam;
    char fname[SEQ_FNAME_MAX_LEN];
    sequence_par *seqp;

    seqp = read_sequence_par("testing_fodder/parameters/sequence.par");
    
    for (cam = 0; cam < 4; cam++) {
        printf("%s", seqp->img_base_name[cam]);
        sprintf(fname, "dumbbell/cam%d_Scene77_", cam + 1);
        fail_unless(strncmp(fname, seqp->img_base_name[cam],
            SEQ_FNAME_MAX_LEN - 1) == 0);
    }
    fail_unless(seqp->first == 497);
    fail_unless(seqp->last == 597);
}
END_TEST

START_TEST(test_read_track_par)
{
    track_par tpar_correct = {
        0.4, 120, 2.0, -2.0, 2.0, -2.0, 2.0, -2.0, 0., 0., 0., 0., 1. 
    };
    
    track_par *tpar;
    tpar = read_track_par("testing_fodder/parameters/track.par");
    
    fail_unless(compare_track_par(tpar, &tpar_correct));
}
END_TEST

START_TEST(test_read_volume_par)
{
    volume_par vpar_correct = {
        {-250., 250.}, {-100., -100.}, {100., 100.}
    };
    
    volume_par *vpar;
    vpar = read_volume_par("testing_fodder/parameters/criteria.par");
    
    fail_unless(compare_volume_par(vpar, &vpar_correct));
}
END_TEST

Suite* fb_suite(void) {
    Suite *s = suite_create ("Parameters handling");

    TCase *tc = tcase_create ("Read sequence parameters");
    tcase_add_test(tc, test_read_sequence_par);
    suite_add_tcase (s, tc);
    
    tc = tcase_create ("Read tracking parameters");
    tcase_add_test(tc, test_read_track_par);
    suite_add_tcase (s, tc);

    tc = tcase_create ("Read illuminated volume parameters");
    tcase_add_test(tc, test_read_volume_par);
    suite_add_tcase (s, tc);

    return s;
}

int main(void) {
    int number_failed;
    Suite *s = fb_suite ();
    SRunner *sr = srunner_create (s);
    srunner_run_all (sr, CK_ENV);
    number_failed = srunner_ntests_failed (sr);
    srunner_free (sr);
    return (number_failed == 0) ? EXIT_SUCCESS : EXIT_FAILURE;
}

