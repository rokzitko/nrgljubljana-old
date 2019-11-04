# Rok Zitko, rok.zitko@ijs.si, May 2014

# CHANGE LOG
# 21.5.2014 - first version

if TOOLS
if HAVE_GSL
bin_PROGRAMS = integ
endif
endif

integ_SOURCES = integ.cc
integ_CXXFLAGS = $(AM_CXXFLAGS) $(GSL_CFLAGS)
integ_LDFLAGS = $(GSL_LIBS)
integ_LDADD = $(GSL_LIBS)

EXTRA_DIST = README
