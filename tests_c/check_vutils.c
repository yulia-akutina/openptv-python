/*  Unit tests for the vector utilities. Uses the Check
    framework: http://check.sourceforge.net/
    
    To run it, type "make check" when in the top C directory, src_c/
    If that doesn't run the tests, use the Check tutorial:
    http://check.sourceforge.net/doc/check_html/check_3.html
*/

#include <check.h>
#include <stdlib.h>

#include "../src_c/vec_utils.h"

START_TEST(test_init_pos3d)
{
    int i;
    pos3d p;
    init_pos3d(p);
    
    for (i = 0; i < 3; i++) fail_unless(p[i] == -999);
}
END_TEST

START_TEST(test_copy_pos3d)
{
    int i;
    pos3d src = {1., 2., 3.}, dst;
    copy_pos3d(dst, src);
    
    for (i = 0; i < 3; i++) fail_unless(dst[i] == src[i]);
}
END_TEST

START_TEST(test_subtract_pos3d)
{
    int i;
    pos3d sub = {1., 2., 3.}, from = {4., 5., 6.}, out;
    subst_pos3d(from, sub, out);
    
    for (i = 0; i < 3; i++) fail_unless(out[i] == 3.);
}
END_TEST

START_TEST(test_diff_norm)
{
    int i;
    pos3d vec1 = {1., 2., 3.}, vec2 = {4., 5., 6.};
    fail_unless(diff_norm_pos3d(vec1, vec2) == sqrt(2)*3);
}
END_TEST

Suite* vu_suite(void) {
    Suite *s = suite_create ("Vector utilities");
    TCase *tc;
    
    tc = tcase_create("Init pos3d");
    tcase_add_test(tc, test_init_pos3d);
    suite_add_tcase (s, tc);
    
    tc = tcase_create("Copy pos3d");
    tcase_add_test(tc, test_copy_pos3d);
    suite_add_tcase (s, tc);

    tc = tcase_create("Subtract pos3d");
    tcase_add_test(tc, test_subtract_pos3d);
    suite_add_tcase (s, tc);

    return s;
}

int main(void) {
    int number_failed;
    Suite *s = vu_suite ();
    SRunner *sr = srunner_create (s);
    srunner_run_all (sr, CK_ENV);
    number_failed = srunner_ntests_failed (sr);
    srunner_free (sr);
    return (number_failed == 0) ? EXIT_SUCCESS : EXIT_FAILURE;
}

