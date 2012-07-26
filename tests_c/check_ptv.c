
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
        
        /* Clear the tracks and leave the testing directory clean: */
        if ((strncmp(next_file->d_name, "added", 5) == 0) || \
            (strncmp(next_file->d_name, "ptv_is", 6) == 0)) {
            remove(res_dir_name);
        }
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
    Exterior correct_ext = {
        105.2632, 102.7458, 403.8822,
        -0.2383291, 0.2442810, 0.0552577, 
        {{0.9688305, -0.0535899, 0.2418587}, 
        {-0.0033422, 0.9734041, 0.2290704},
        {-0.2477021, -0.2227387, 0.9428845}}};
    
    Interior I;
    Interior correct_int = {-2.4742, 3.2567, 100.0000};
    
    Glass G;
    Glass correct_glass = {0.0001, 0.00001, 150.0};
    
    ap_52 addp;
    ap_52 correct_addp = {0., 0., 0., 0., 0., 1., 0.};
    
    Calibration cal;
    Calibration correct_cal = {correct_ext, correct_int, correct_glass, 
        correct_addp};
    
    char ori_file[] = "testing_fodder/cal/cam1.tif.ori";
    char add_file[] = "testing_fodder/cal/cam1.tif.addpar";
    
    fail_unless(read_ori(&Ex, &I, &G, ori_file, &addp, add_file, NULL));
    cal.ext_par = Ex;
    cal.int_par = I;
    cal.glass_par = G;
    cal.added_par = addp;
    
    fail_unless(compare_calib(&cal, &correct_cal));
}
END_TEST

Suite* ptv_suite(void) {
    Suite *s = suite_create ("PTV");

    TCase *tc_tw = tcase_create ("Tracking window");
    tcase_add_test(tc_tw, test_allocate_tracking_structs);
    suite_add_tcase (s, tc_tw);

    TCase *tc_track = tcase_create ("Tracking");
    tcase_add_test(tc_track, test_tracking);
    suite_add_tcase (s, tc_track);

    TCase *tc_rori = tcase_create ("Read orientation file");
    tcase_add_test(tc_rori, test_read_ori);
    suite_add_tcase (s, tc_rori);

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

