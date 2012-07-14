## Process this file with automake to produce Makefile.in

lib_LTLIBRARIES = libptv.la
libptv_la_SOURCES = ptv.c tools.c ttools.c imgcoord.c multimed.c ray_tracing.c\
    lsqadj.c trafo.c ptv.h
libptv_la_LIBADD = -lm

bin_PROGRAMS = main
main_SOURCES = main.c
main_LDADD = libptv.la
