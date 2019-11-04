# Rok Zitko, rok.zitko@ijs.si, May 2014

# CHANGE LOG
# 22.5.2014 - first version

if TOOLS
if HAVE_GSL
bin_PROGRAMS = resample
endif
endif

resample_SOURCES = resample.cc
resample_CXXFLAGS = $(AM_CXXFLAGS) $(GSL_CFLAGS)
resample_LDFLAGS = $(GSL_LIBS)
resample_LDADD = $(GSL_LIBS)

EXTRA_DIST = README
