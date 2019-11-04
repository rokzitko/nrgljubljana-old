# Rok Zitko, rok.zitko@ijs.si, Dec 2009

# CHANGE LOG
# 1.4.2011 - compilation of kk is conditional

if TOOLS
if HAVE_GSL
bin_PROGRAMS = kk
endif
endif

kk_SOURCES = kk.cc
kk_CXXFLAGS = $(AM_CXXFLAGS) $(GSL_CFLAGS)
kk_LDFLAGS = $(GSL_LIBS)
kk_LDADD = $(GSL_LIBS)

EXTRA_DIST = README
