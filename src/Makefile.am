# automake script for the NRG Ljubljana C++ code

bin_PROGRAMS = nrg nrgcmpl

# Make sure src/boost comes before system include. We want
# to use bundled boost ublas headers, otherwise all hell
# breaks loose.

AM_CXXFLAGS = -I$(top_srcdir)/coefnew \
  -I$(top_srcdir)/src -I$(top_srcdir)/src/boost \
  $(BOOST_CPPFLAGS) $(PTHREAD_CFLAGS)

AM_LDFLAGS  = $(BOOST_LDFLAGS)

LDADD = $(LAPACK_LIBS) $(BLAS_LIBS) $(MKL_LIBS) $(MISC_LIBS) $(LIBS) $(PTHREAD_LIBS) \
  $(BOOST_SERIALIZATION_LIB) $(BOOST_MPI_LIB)

# -DNDEBUG removes the debugging code of boost::ublas
AM_CXXFLAGS += -DNDEBUG

nrg_SOURCES = nrg.cc nrg.h $(nrg_coefs)
nrg_CXXFLAGS = -DNRG_REAL $(AM_CXXFLAGS)

nrgcmpl_SOURCES = nrg.cc nrg.h $(nrg_coefs)
nrgcmpl_CXXFLAGS = -DNRG_COMPLEX $(AM_CXXFLAGS)

noinst_HEADERS = nrg-lr-common.cc portabil.h spectral.h debug.h \
misc.h numerics.h bins.h dmnrg.h tridiag.h diag.h \
symmetry.cc matrix.cc param.cc outfield.cc invar.cc \
read-input.cc recalc.cc spec.cc splitting.cc \
spec_FT.cc spec_DMNRG.cc spec_FDM.cc spec_CFS.cc spectrumrealfreq.cc \
sym-ANYJ.cc sym-DBLISOSZ.cc sym-DBLSU2.cc sym-ISO.cc sym-ISOLR.cc \
sym-ISOSZ.cc sym-ISOSZLR.cc sym-NONE.cc sym-QS.cc sym-QSC3.cc \
sym-QSLR.cc sym-QST.cc sym-QSZ.cc sym-QSZLR.cc sym-SL.cc sym-SL3.cc \
sym-SPSU2.cc sym-SPSU2C3.cc sym-SPSU2LR.cc sym-SPSU2T.cc sym-SPU1.cc \
sym-SPU1LR.cc sym-SU2.cc sym-U1.cc sym-QSTZ.cc sym-QSZTZ.cc sym-QJ.cc \
sym-P.cc sym-PP.cc \
nrg-recalc-QS.cc nrg-recalc-QSLR.cc nrg-recalc-QSZ.cc nrg-recalc-QSZLR.cc \
nrg-recalc-ISO.cc nrg-recalc-ISOLR.cc nrg-recalc-U1.cc\
nrg-recalc-NONE.cc nrg-recalc-ANYJ.cc \
nrg-recalc-ANYJ.cc nrg-recalc-QSLR.cc nrg-recalc-SL.cc nrg-recalc-ANYJ.cc.m4 \
nrg-recalc-ISO.cc.m4 nrg-recalc-ISOLR.cc.m4 nrg-recalc-NONE.cc.m4 \
nrg-recalc-QS.cc.m4 nrg-recalc-QSLR.cc.m4 nrg-recalc-QSZ.cc.m4 \
nrg-recalc-QSZLR.cc.m4 nrg-recalc-SL.cc.m4 nrg-recalc-U1.cc.m4 \
nrg-recalc-SPSU2.cc nrg-recalc-SPSU2.cc.m4 nrg-recalc-SPU1.cc \
nrg-recalc-SPU1.cc.m4 nrg-recalc-ISO2.cc \
nrg-recalc-ISO2.cc.m4 nrg-recalc-ISO2LR.cc nrg-recalc-ISO2LR.cc.m4 \
nrg-recalc-ISOSZ.cc nrg-recalc-ISOSZLR.cc nrg-recalc-SU2.cc.m4 \
nrg-recalc-SU2.cc nrg-recalc-DBLSU2.cc nrg-recalc-DBLSU2.cc.m4 \
nrg-recalc-DBLISOSZ.cc nrg-recalc-SPU1LR.cc \
nrg-recalc-SL3.cc.m4 nrg-recalc-SL3.cc \
nrg-recalc-ISOSZ.cc.m4 nrg-recalc-ISOSZLR.cc.m4 \
nrg-recalc-SPU1LR.cc.m4 nrg-recalc-DBLISOSZ.cc.m4 \
nrg-recalc-SPSU2LR.cc.m4 nrg-recalc-SPSU2LR.cc \
nrg-recalc-QST.cc.m4 nrg-recalc-QST.cc \
nrg-recalc-SPSU2T.cc.m4 nrg-recalc-SPSU2T.cc \
nrg-recalc-QSC3.cc.m4 nrg-recalc-QSC3.cc \
nrg-recalc-SPSU2C3.cc.m4 nrg-recalc-SPSU2C3.cc \
nrg-recalc-QSTZ.cc.m4 nrg-recalc-QSTZ.cc \
nrg-recalc-QSZTZ.cc.m4 nrg-recalc-QSZTZ.cc \
nrg-recalc-QJ.cc.m4 nrg-recalc-QJ.cc \
nrg-recalc-P.cc.m4 nrg-recalc-P.cc \
nrg-recalc-PP.cc.m4 nrg-recalc-PP.cc \
matsubara.h matsubara2.h

EXTRA_DIST = rational.hpp m4-update recalc-macros.m4 cblas.h boost/
