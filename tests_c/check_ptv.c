
/*  Unit-test suit for the C core of PyPTV. Uses the Check framework:
    http://check.sourceforge.net/
    
    To run it, type "make check" when in the top C directory, src_c/
*/

#include <check.h>
#include <stdlib.h>

#include <stdio.h>
#include <unistd.h>
#include <sys/types.h>
#include <dirent.h>

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

/* Check that target t1 is equal to target t2, i.e. all their fields are equal.
*/
int compare_targets(target *t1, target *t2) {
    return (\
        (t1->pnr == t2->pnr) && (t1->x == t2->x) && (t1->y == t2->y) && \
        (t1->n == t2->n) && (t1->nx == t2->nx) && (t1->ny == t2->ny) && \
        (t1->sumg == t2->sumg) && (t1->tnr == t2->tnr));
}

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

int compare_corres(corres *c1, corres *c2) {
    return ((c1->nr == c2->nr) && (c1->p[0] == c2->p[0]) && \
        (c1->p[1] == c2->p[1]) && (c1->p[2] == c2->p[2]) && \
        (c1->p[3] == c2->p[3]));
}

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

/* Test a full tracking run with trackcor_*, from init to completion.
 * This is done by running on an initial set of rt_is files, then comparing the
 * result to a sample run.
 */
START_TEST(test_tracking)
{
    int step = 0;
    DIR *res_dir = NULL;
    struct dirent *next_file;
    FILE *res_file, *res_sample;
    char res_dir_name[128] = "res/", samp_dir_name[128] = "sample_res/";
    char **line_res, **line_samp;
    size_t buf_size = 128;
    int line_len = 0;
    
    /* For now, the code relies on a rigid experiment directory structure and
       expects to find there the input files and parameters. It also throws the
       output in the same tree. Until we fix that, testing_fodder/ mimics the 
       necessary structure. */
    fail_unless(!chdir("testing_fodder/"));
    
    trackcorr_c_init();
    for (step = 497; step < 597; step++) {
        trackcorr_c_loop(step, lmax_track, ymin_track, ymax_track, 0);
    }
    trackcorr_c_finish(597);
    
    /* After tracking all outputs are in res/, and compared against 
       sample_res/
    */
    res_dir = opendir("res/");
    fail_if(res_dir == NULL);
    
    line_res = (char**) malloc(buf_size);
    line_samp = (char**) malloc(buf_size);
    
    while (next_file = readdir(res_dir)) {
        strcpy(res_dir_name + 4, next_file->d_name);
        fail_if((res_file = fopen(res_dir_name, "r")) == NULL);
        
        strcpy(samp_dir_name + 11, next_file->d_name);
        fail_if((res_sample = fopen(res_dir_name, "r")) == NULL);
        
        while ((line_len = getline(line_res, &buf_size, res_file)) != -1) {
            fail_unless(line_len = getline(line_samp, &buf_size, res_sample));
            fail_if(strncmp(*line_res, *line_samp, line_len));
        }
        fclose(res_file);
        fclose(res_sample);
    }
    closedir(res_dir);
}
END_TEST

/* Regression test for reading orientation files. Just reads a sample file and
   makes sure that nothing crashes and the orientation structures are filled
   out correctly.
*/
START_TEST(test_read_ori)
{
    Exterior Ex;
    Interior I;
    Glass    G;
    ap_52    addp;
    char     ori_file[] = "cal/cam1.tif.ori", add_file[] = , add_fallback[] = ;
    
    fail_unless(read_ori(Ex, I, G, ori_file, addp, add_file, add_fallback));
}
END_TEST

Suite* ptv_suite(void) {
    Suite *s = suite_create ("PTV");

    TCase *tc_tw = tcase_create ("Tracking window");
    tcase_add_test(tc_tw, test_allocate_tracking_structs);
    suite_add_tcase (s, tc_tw);
    
    TCase *tc_trt = tcase_create ("Read targets");
    tcase_add_test(tc_trt, test_read_targets);
    suite_add_tcase (s, tc_trt);

    TCase *tc_trpf = tcase_create ("Read path frame");
    tcase_add_test(tc_trpf, test_read_path_frame);
    suite_add_tcase (s, tc_trpf);

    TCase *tc_track = tcase_create ("Tracking");
    tcase_add_test(tc_track, test_tracking);
    suite_add_tcase (s, tc_track);
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

