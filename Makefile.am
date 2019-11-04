ACLOCAL_AMFLAGS = -I m4

# Include major version number
pkgdatadir = ${datarootdir}/$(PACKAGE)-$(PACKAGE_VERSION)

# NOTE: "Such scripts for which a build rule has been supplied need to be
# deleted explicitly using CLEANFILES (see Clean), and their sources have to
# be distributed, usually with EXTRA_DIST (see Basics of Distribution)."

nrginit: $(srcdir)/nrginit.sh
	rm -f nrginit
	echo "#!"$(BASH) >nrginit
	m4 -DPKGDATADIR=$(pkgdatadir) $(srcdir)/nrginit.sh >>nrginit
	chmod ugo+x nrginit

nrginit.m: $(srcdir)/nrginit.m.m4
	m4 -DPKGDATADIR=$(pkgdatadir) -DVERSION=$(PACKAGE_VERSION) \
	    $(srcdir)/nrginit.m.m4 >nrginit.m

# The nrgrun* scripts are mostly for backwards compatibility. It's now
# recommended to call the executables directly.

nrgrun: nrgrun.sh
	m4 -DVERSION=$(PACKAGE_VERSION) nrgrun.sh >nrgrun
	chmod a+x nrgrun

bin_SCRIPTS = nrginit nrgrun
CLEANFILES = $(bin_SCRIPTS) $(check_SCRIPTS) nrginit.m

EXTRA_DIST = reconf configure nrginit.sh nrginit.m.m4 README.Mac nrgrun.sh

dist_pkgdata_DATA = sneg.m initial.m models.m \
                    initialparse.m operators.m dmft.m
pkgdata_DATA = nrginit.m

SUBDIRS = m4 coefnew src scripts tools
