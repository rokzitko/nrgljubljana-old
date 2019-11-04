# Decoupling tools
# Rok Zitko, rok.zitko@ijs.si, Mar 2011

AM_CXXFLAGS = -I$(top_srcdir)/coefnew \
  -I$(top_srcdir)/src -I$(top_srcdir)/src/boost \
   $(BOOST_CPPFLAGS) $(PTHREAD_CFLAGS)
    
AM_LDFLAGS  = $(BOOST_LDFLAGS)
    
diag_LDADD = $(MKL_LIBS) $(MISC_LIBS)

if TOOLS
bin_PROGRAMS = diag unitary matrix
endif

diag_SOURCES = diag.cc

unitary_SOURCES = unitary.cc

AM_YFLAGS=-d

matrix_SOURCES = parser.yy matrix.ll matrix.h

EXTRA_DIST = README
