## Process this file with automake to produce Makefile.in

check_LTLIBRARIES = libptv.la libpar.la libvecutils.la

libptv_la_SOURCES = tools.c ttools.c imgcoord.c multimed.c ray_tracing.c\
    lsqadj.c trafo.c ptv.h track.c pointpos.c intersect.c orientation.c \
    rotation.c epi.c calibration.c \
    vec_utils.c vec_utils.h jw_ptv.c segmentation.c peakfitting.c \
    correspondences.c sortgrid.c image_processing.c parameters.c parameters.h \
    tracking_run.c tracking_run.h corresp.c btree.c
libptv_la_LIBADD = -lm -loptv

libpar_la_SOURCES = parameters.c parameters.h

libvecutils_la_SOURCES = vec_utils.c vec_utils.h
libvecutils_la_LIBADD = -lm

check_PROGRAMS = main
main_SOURCES = main.c
main_LDADD = libptv.la
