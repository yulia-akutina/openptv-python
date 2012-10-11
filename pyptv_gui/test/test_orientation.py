# Regression tests for the orienbtation.

import unittest
import os, shutil
import numpy as np

from ptv1 import py_start_proc_c, py_init_proc_c, py_prepare_eval

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
        shutil.copyfile("parameters/sequence_scene.par", "parameters/sequence.par")
        os.remove("parameters/sequence_scene.par")
        os.chdir("../")
    
    def test_prepare_eval(self):
        shutil.copytree("db_targ/", "scene83_event1/")
        shutil.copytree("db_res/", "res/")
        shutil.copyfile("parameters/sequence.par", "parameters/sequence_scene.par")
        shutil.copyfile("parameters/sequence_db.par", "parameters/sequence.par")
        
        py_init_proc_c()
        py_start_proc_c()
        
        pix, crd = py_prepare_eval(4)
        
        #np.savetxt('db_res/pix_regress.dat', pix.reshape(-1, 2))
        #np.savetxt('db_res/crd_regress.dat', crd.reshape(-1, 3))
        np.testing.assert_array_equal(pix.reshape(-1, 2),
            np.loadtxt('db_res/pix_regress.dat'))
        np.testing.assert_array_equal(crd.reshape(-1, 3),
            np.loadtxt('db_res/crd_regress.dat'))
        
if __name__ == '__main__':
    unittest.main()

