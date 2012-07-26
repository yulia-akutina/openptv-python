/*  Unit tests for functions related to the frame buffer. Uses the Check
    framework: http://check.sourceforge.net/
    
    To run it, type "make check" when in the top C directory, src_c/
    If that doesn't run the tests, use the Check tutorial:
    http://check.sourceforge.net/doc/check_html/check_3.html
*/

#include <check.h>
#include <stdlib.h>
#include <stdio.h>

#include "../src_c/tracking_frame_buf.h"

START_TEST(test_read_targets)
{
    target tbuf[2]; /* Two targets in the sample target file */
    target t1 = {0, 1127.0000, 796.0000, 13320, 111, 120, 828903, 1};
    target t2 = {1, 796.0000, 809.0000, 13108, 113, 116, 658928, 0};
    
    char *file_base = "testing_fodder/sample_";
    int frame_num = 42;
    int targets_read = 0;
    
    targets_read = read_targets(tbuf, file_base, frame_num);
    fail_unless(targets_read == 2);
    fail_unless(compare_targets(tbuf, &t1));
    fail_unless(compare_targets(tbuf + 1, &t2));
}
END_TEST

START_TEST(test_write_targets)
{
    /* Write and read back two targets, make sure they're the same.
    Assumes that read_targets is ok, which is tested separately. */
    
    target tbuf[2]; /* Two targets in the sample target file */
    target t1 = {0, 1127.0000, 796.0000, 13320, 111, 120, 828903, 1};
    target t2 = {1, 796.0000, 809.0000, 13108, 113, 116, 658928, 0};
    
    char *file_base = "testing_fodder/test_";
    int frame_num = 42;
    int num_targets = 2;
    
    tbuf[0] = t1; tbuf[1] = t2;
    fail_unless(write_targets(tbuf, num_targets, file_base, frame_num));
    
    num_targets = read_targets(tbuf, file_base, frame_num);
    fail_unless(num_targets == 2);
    fail_unless(compare_targets(tbuf, &t1));
    fail_unless(compare_targets(tbuf + 1, &t2));
    
    // Leave the test directory clean:
    remove("testing_fodder/test_0042_targets");
}
END_TEST

START_TEST(test_read_path_frame)
{
    corres cor_buf[80];
    P path_buf[80];
    int alt_link;
    
    /* Correct values for particle 3 */
    P path_correct = {
        .x = {45.219, -20.269, 25.946},
        .prev = -1,
        .next = -2,
        .prio = 4, 
        .finaldecis = 1000000.0,
        .inlist = 0.
    };
    for (alt_link = 0; alt_link < POSI; alt_link++) {
        path_correct.decis[alt_link] = 0.0;
        path_correct.linkdecis[alt_link] = -999;
    }
    corres c_correct = { 3, {96, 66, 26, 26} };

    char *file_base = "testing_fodder/rt_is";
    int frame_num = 818;
    int targets_read = 0;

    targets_read = read_path_frame(cor_buf, path_buf, file_base, frame_num);
    fail_unless(targets_read == 80);
    fail_unless(compare_corres(cor_buf + 2, &c_correct), \
        "Got corres: %d, [%d %d %d %d]", cor_buf[2].nr, \
        cor_buf[2].p[0], cor_buf[2].p[1], cor_buf[2].p[2], cor_buf[2].p[3]);
    fail_unless(compare_path_info(path_buf + 2, &path_correct));
}
END_TEST

Suite* fb_suite(void) {
    Suite *s = suite_create ("Frame Buffer");

    TCase *tc_trt = tcase_create ("Read targets");
    tcase_add_test(tc_trt, test_read_targets);
    suite_add_tcase (s, tc_trt);

    TCase *tc_twt = tcase_create ("Write targets");
    tcase_add_test(tc_twt, test_write_targets);
    suite_add_tcase (s, tc_twt);

    TCase *tc_trpf = tcase_create ("Read path frame");
    tcase_add_test(tc_trpf, test_read_path_frame);
    suite_add_tcase (s, tc_trpf);

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

