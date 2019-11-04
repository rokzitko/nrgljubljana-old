# Hilbert transformation tool
# Rok Zitko, rok.zitko@ijs.si, Dec 2012

if TOOLS
if HAVE_GSL
bin_PROGRAMS = hilb
endif
endif

hilb_SOURCES = hilb.cc
hilb_CXXFLAGS = $(AM_CXXFLAGS) $(GSL_CFLAGS)
hilb_LDFLAGS = $(GSL_LIBS)
hilb_LDADD = $(GSL_LIBS)

EXTRA_DIST = README

