## Process this file with automake to produce Makefile.in

lib_LTLIBRARIES = libptv.la libfb.la libpar.la libvecutils.la

libptv_la_SOURCES = tools.c ttools.c imgcoord.c multimed.c ray_tracing.c\
    lsqadj.c trafo.c ptv.h track.c pointpos.c intersect.c orientation.c \
    rotation.c epi.c calibration.c tracking_frame_buf.c tracking_frame_buf.h \
    vec_utils.c vec_utils.h jw_ptv.c segmentation.c peakfitting.c \
    correspondences.c sortgrid.c image_processing.c parameters.c parameters.h \
    tracking_run.c tracking_run.h corresp.c btree.c btree.h
libptv_la_LIBADD = -lm

libfb_la_SOURCES = tracking_frame_buf.c tracking_frame_buf.h

libpar_la_SOURCES = parameters.c parameters.h

libvecutils_la_SOURCES = vec_utils.c vec_utils.h
libvecutils_la_LIBADD = -lm

bin_PROGRAMS = main
main_SOURCES = main.c
main_LDADD = libptv.la
