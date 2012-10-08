# Regression tests for the orienbtation.

import unittest
import os, shutil
import numpy as np

from ptv1 import py_start_proc_c, py_init_proc_c

class TestOrient(unittest.TestCase):
    def setUp(self):
        os.chdir("testing_fodder/")

        if os.path.exists("res/"):
            shutil.rmtree("res/")
        if os.path.exists("scene83_event1/"):
            shutil.rmtree("scene83_event1/")
        
    def tearDown(self):
        shutil.rmtree("res/")
        shutil.rmtree("scene83_event1/")
        os.chdir("../")

    def test_prepare_eval(self):
        shutil.copytree("db_targ/", "scene83_event1/")
        shutil.copytree("db_res/", "res/")
        
        py_init_proc_c()
        py_start_proc_c()
        
        pix, crd = py_prepare_eval(4)
        
        np.savetxt('pix_regress.dat', pix)
        np.savetxt('crd_regress.dat', crd)
        #np.testing.assert_array_equal(pix, np.loadtxt('pix_regress.dat'))
        #np.testing.assert_array_equal(crd, np.loadtxt('crd_regress.dat'))
        
if __name__ == '__main__':
    unittest.main()

