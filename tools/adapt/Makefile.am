# automake script for the 'adapt' tools for adaptable discretization grid
# Rok Zitko, rok.zitko@ijs.si

# CHANGE LOG
# 1.4.2011 - compilaiton of nrgchain is conditional

if TOOLS
bin_PROGRAMS = adapt
endif

if TOOLS
# nrgchain requires libgmp
if HAVE_GMP
bin_PROGRAMS += nrgchain
endif
endif

#LDADD = $(LIBS) $(PTHREAD_LIBS)

adapt_SOURCES = adapt.cc

nrgchain_SOURCES = nrgchain.cc

noinst_HEADERS = calc.h  io.h  lambda.h  linint.h  load.h  parser.h

EXTRA_DIST = README param Delta.dat FSOL.dat.ref FSOLNEG.dat.ref \
             GSOL.dat.ref GSOLNEG.dat.ref xi.dat.ref zeta.dat.ref

quicktest:
	./adapt P
	./adapt N
	./nrgchain
	diff FSOL.dat FSOL.dat.ref
	diff FSOLNEG.dat FSOLNEG.dat.ref
	diff GSOL.dat GSOL.dat.ref
	diff GSOLNEG.dat GSOLNEG.dat.ref
	diff xi.dat xi.dat.ref
	diff zeta.dat zeta.dat.ref
