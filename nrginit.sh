# Start Mathematica in batch mode and run the "NRG Ljubljana"
# initialisation script
# Rok Zitko, rok.zitko@ijs.si, Dec 2009, Nov 2013

# The path is hard-coded during the installation of "NRG Ljubljana".
NRG_INITSCRIPT=PKGDATADIR/nrginit.m

math -batchinput -batchoutput <${NRG_INITSCRIPT}
