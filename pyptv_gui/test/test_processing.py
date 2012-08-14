# Regression tests for the post-processing tests: sequencing, tracking, tracking 
# back, etc.

import unittest
import os, shutil, re

from scipy.misc import imread

from ptv1 import py_start_proc_c, py_init_proc_c
from ptv1 import py_sequence_init, py_sequence_loop, py_set_img

def compare_directories(dir1, dir2, mask=None):
    """
    For all files in dir1, check that the file exists in dir2 and that the
    files are equal line by line.
    
    Arguments:
    dir1, dir2 - path to the directories for comparison.
    mask - a compiled regular wxpression. Only file names matching it will be
        compared, others will be skipped. If None (default), all files are
        compared.
    
    Returns:
    True if equal, False otherwise.
    """
    for fname in os.listdir(dir1):
        if (mask is not None) and (mask.search(fname) is None):
            continue
        
        fout = open(dir1 + fname, "r")
        fcmp = open(dir2 + fname, "r")
        
        while True:
            lout = fout.readline()
            lcmp = fcmp.readline()
            if lout != lcmp:
                return False
            if len(lout) == 0: break
        
        fout.close()
        fcmp.close()
    
    return True

class TestProcessing(unittest.TestCase):
    def setUp(self):
        os.chdir("testing_fodder/")

        if os.path.exists("res/"):
            shutil.rmtree("res/")
        if os.path.exists("scene83_event1/"):
            shutil.rmtree("scene83_event1/")
        
        shutil.copytree("clean_scene/", "scene83_event1/")
        
    def tearDown(self):
        shutil.rmtree("res/")
        shutil.rmtree("scene83_event1/")
    
    def test_sequencing(self):
        """Sequencing reproduces sample results"""
        os.mkdir("res/")
        
        py_init_proc_c()
        py_start_proc_c()
        py_sequence_init(0)
        
        for frame in xrange(497, 598):
            for cam in xrange(4):
                img = imread(
                    "scene83_event1/cam%d_Scene83_%04d" % (cam + 1, frame))
                py_set_img(img, cam)
                py_sequence_loop(0, frame)
        
        self.failUnless(compare_directories("res/", "after_sequencing/"))
        self.failUnless(compare_directories(
            "scene83_event1/", "after_sequencing_targets/",
            mask=re.compile("_targets$")))

