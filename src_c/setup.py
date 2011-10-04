# -*- coding: utf-8 -*-
# sqadj.c pointpos.c multimed.c imgcoord.c intersect.c ray_tracing.c trafo.c
from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext

import numpy as np

setup(
    name="ptv1",
    cmdclass = {'build_ext': build_ext},
    ext_modules = [Extension("ptv1", ["ptv1.pyx", "segmentation.c", "tools.c","image_processing.c", "trafo.c", "jw_ptv.c", "peakfitting.c", "rotation.c", "correspondences.c", "epi.c", "multimed.c", "ray_tracing.c","imgcoord.c","lsqadj.c", "orientation.c","sortgrid.c", "pointpos.c","intersect.c","track.c","ttools.c","ptv.c","draw.c","mousefunction.c"],
                             include_dirs = [np.get_include(),'.'],
                             extra_compile_args=['-O3'])],
    py_modules = ['ptv1',],
)

