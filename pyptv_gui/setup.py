# -*- coding: utf-8 -*-
# sqadj.c pointpos.c multimed.c imgcoord.c intersect.c ray_tracing.c trafo.c
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
        "mousefunction.c", "tracking_frame_buf.c", "vec_utils.c"],
        include_dirs = [np.get_include(),'.'],
        extra_compile_args=['-O3', '-fno-common']),
    Extension("tracking_framebuf", ["tracking_framebuf.pyx",
        "tracking_frame_buf.c"]) ]

setup(
    name="ptv1",
    cmdclass = {'build_ext': build_ext},
    ext_modules = ext_mods,
    py_modules = ['ptv1',],
)


