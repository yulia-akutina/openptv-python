# A test for tracking backwards. Tracking forward is tested at the C level
# right now, although in the future we'd want to test it more, and this would
# be the place for it.

import unittest
import os, shutil

from ptv1 import py_start_proc_c, py_trackback_c, py_init_proc_c

class TestTrackback(unittest.TestCase):
    def setUp(self):
        os.chdir("testing_fodder/")
        
        if os.path.exists("res/"):
            shutil.rmtree("res/")
        shutil.copytree("sample_res/", "res/")
    
    def tearDown(self):
        shutil.rmtree("res/")
    
    def test_backwards(self):
        """Backward tracking reproduces sample results"""
        py_init_proc_c()
        py_start_proc_c()
        py_trackback_c()
        
        for fname in os.listdir("res/"):
            print fname
            fout = open("res/" + fname, "r")
            fcmp = open("backward_res/" + fname, "r")
            
            while True:
                lout = fout.readline()
                lcmp = fcmp.readline()
                self.failUnlessEqual(lout, lcmp, "File %s: %s != %s" % (fname, lout.strip(), lcmp.strip()))
                if len(lout) == 0: break
            
            fout.close()
            fcmp.close()

