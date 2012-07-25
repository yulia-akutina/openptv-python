## Process this file with automake to produce Makefile.in

lib_LTLIBRARIES = libptv.la libfb.la
libptv_la_SOURCES = ptv.c tools.c ttools.c imgcoord.c multimed.c ray_tracing.c\
    lsqadj.c trafo.c ptv.h track.c pointpos.c intersect.c orientation.c \
    rotation.c epi.c calibration.c tracking_frame_buf.c
libptv_la_LIBADD = -lm

libfb_la_SOURCES = tracking_frame_buf.c tracking_frame_buf.h

bin_PROGRAMS = main
main_SOURCES = main.c
main_LDADD = libptv.la
