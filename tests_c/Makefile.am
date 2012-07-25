TESTS = check_ptv check_fb
check_PROGRAMS = check_ptv check_fb

check_ptv_SOURCES = check_ptv.c $(top_builddir)/src_c/ptv.h
check_ptv_CFLAGS = @CHECK_CFLAGS@
check_ptv_LDADD = $(top_builddir)/src_c/libptv.la @CHECK_LIBS@ -lm

check_fb_SOURCES = check_fb.c $(top_builddir)/src_c/tracking_frame_buf.h
check_fb_CFLAGS = @CHECK_CFLAGS@
check_fb_LDADD = $(top_builddir)/src_c/libfb.la @CHECK_LIBS@
