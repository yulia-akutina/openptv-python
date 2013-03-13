# -*- coding: utf-8 -*-
from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext

import numpy as np
import os

import sys
if sys.platform.startswith(("linux","darwin")):
    output_name = 'ptv1.so'
else:
    output_name = 'ptv1.pyd'

origdir = os.getcwd()
os.chdir('../src_c')

ext_mods = [
    Extension("ptv1", ["ptv1.pyx", "segmentation.c", "tools.c",
        "image_processing.c", "trafo.c", "jw_ptv.c", "peakfitting.c", 
        "rotation.c", "correspondences.c", "epi.c", "multimed.c", 
        "ray_tracing.c", "imgcoord.c", "lsqadj.c", "orientation.c","sortgrid.c",
        "pointpos.c", "intersect.c", "track.c", "ttools.c", "draw.c",
        "mousefunction.c", "vec_utils.c", "parameters.c", "tracking_run.c"],
        include_dirs=[np.get_include(),'.'], libraries=['optv'],
        extra_compile_args=['-O3', '-fno-common']),
    Extension("tracking_framebuf", ["tracking_framebuf.pyx"], libraries=['optv']),
    Extension("tracking_run_py", ["tracking_run_py.pyx", "tracking_run_py.pxd", 
        "tracking_run.c", "parameters.c"], libraries=['optv']) ]

setup(
    name="ptv1",
    cmdclass = {'build_ext': build_ext},
    ext_modules = ext_mods,
    py_modules = ['ptv1',],
)


