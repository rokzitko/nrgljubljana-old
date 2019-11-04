# autoconf script for "NRG Ljubljana"
# Rok Zitko, rok.zitko@ijs.si

AC_INIT(NRG Ljubljana, 2.4.4.25, Rok Zitko <rok.zitko@ijs.si>, nrgljubljana)
AC_CONFIG_MACRO_DIR(m4)
AM_CONFIG_HEADER(config.h)
AM_INIT_AUTOMAKE

# Check for programs
AC_PROG_CC
AC_PROG_CPP
AC_AIX
AC_ISC_POSIX
AC_MINIX
AC_HEADER_STDC

AC_PROG_CXX
AC_PROG_CXXCPP

AC_PROG_INSTALL
AC_PATH_PROGS(BASH, bash)
AC_PATH_PROGS(PERL, perl perl5.004 perl5.003 perl5.002 perl5.001 perl5)

AM_PATH_GSL(1.0, [], [])

# Check for libraries
AX_PTHREAD

AC_CHECK_HEADERS([execinfo.h], [], [], [])

# Did we ask for mkl
AC_ARG_WITH([mkl],
  # Help string
  AC_HELP_STRING([--with-mkl], [use MKL]),
  # action if given
  [
    AC_DEFINE([MKL], [yes], [MKL Libraries])
    MKL_LIBS="$MKL_LIBS"
    AC_SUBST(MKL_LIBS)
  ], [])

AC_ARG_WITH([misc],
  # Help string
  AC_HELP_STRING([--with-misc], [Additional libraries]),
  # action if given
  [
    AC_DEFINE([MISC], [yes], [Misc])
    MISC_LIBS="$MISC_LIBS"
    AC_SUBST(MISC_LIBS)
  ], [])

# Compile tools?
AC_ARG_WITH([tools],
  # Help string
  AC_HELP_STRING([--with-tools], [compile and install tools]),
  # action if given
  [withtools=yes], [])

AM_CONDITIONAL(TOOLS, test x$withtools = xyes)

AX_BOOST_BASE(1.34)
AX_BOOST_SERIALIZATION(1.34)
AX_BOOST_MPI(1.41)
#AX_BLAS([], [AC_WARN(Cannot find BLAS)])
#AX_LAPACK([], [AC_WARN(Cannot find LAPACK)])

AC_CHECK_LIB(gmp, __gmpz_init,
[
 AC_DEFINE([HAVE_LIBGMP],[yes],[GNU Multiple precision libraries])
 LIBS="$LIBS -lgmp"
 havelibgmp=true
],
[
 AC_WARN([GMP not found. Tridiagonalisation code disabled.])
 havelibgmp=false
])

AM_CONDITIONAL(HAVE_GMP, test x$havelibgmp = xtrue)

AC_CHECK_LIB(rt, clock_gettime, ,
[AC_WARN([librt not found. Timing code disabled.])])

AC_CHECK_FUNCS(sqr)

# For the decoupling tools
AC_PROG_YACC
AC_PROG_LEX

AC_ARG_WITH(inspect,
  [ AC_HELP_STRING([--with-inspect], [Inspect test suite errors]) ],
  [ AC_DEFINE(INSPECT_ERRORS, 1, "Inspect test suite errors")],
  [ AC_DEFINE(INSPECT_ERRORS, 0, "Inspect test suite errors")])

AC_CONFIG_FILES(Makefile m4/Makefile \
coefnew/Makefile coefnew/qs/Makefile coefnew/qsz/Makefile \
coefnew/iso/Makefile coefnew/qslr/Makefile coefnew/qszlr/Makefile \
coefnew/isolr/Makefile coefnew/none/Makefile coefnew/u1/Makefile \
coefnew/iso2/Makefile coefnew/iso2lr/Makefile \
coefnew/spsu2/Makefile coefnew/spu1/Makefile coefnew/spu1lr/Makefile \
coefnew/sl/Makefile coefnew/anyj/Makefile \
coefnew/isosz/Makefile coefnew/isoszlr/Makefile \
coefnew/su2/Makefile coefnew/dblsu2/Makefile \
coefnew/dblisosz/Makefile coefnew/sl3/Makefile \
coefnew/qst/Makefile coefnew/spsu2lr/Makefile \
coefnew/spsu2t/Makefile \
coefnew/qsc3/Makefile \
coefnew/spsu2c3/Makefile \
coefnew/qstz/Makefile \
coefnew/qsztz/Makefile \
coefnew/qj/Makefile \
coefnew/p/Makefile \
coefnew/pp/Makefile \
src/Makefile \
scripts/Makefile \
tools/Makefile \
tools/kk-GSL/Makefile \
tools/adapt/Makefile \
tools/bw/Makefile \
tools/broaden/Makefile \
tools/mats/Makefile \
tools/intavg/Makefile \
tools/tdavg/Makefile \
tools/dec/Makefile \
tools/hilb/Makefile \
tools/integrate/Makefile \
tools/resample/Makefile )
AC_OUTPUT
