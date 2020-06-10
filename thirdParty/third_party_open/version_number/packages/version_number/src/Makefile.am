#-------------------------------------------------------------------------------
#   Makefile for version_number
#
#   Adri.Mourits@deltares.nl
#
#   Copyright (C)  Stichting Deltares, 2011-2013.
#-------------------------------------------------------------------------------

include $(top_srcdir)/common.am
noinst_nodist_PROGRAMS = version_number
BUILT_SOURCES = version_number

#-------------------------------------------------------------------------------
#Includes

#-------------------------------------------------------------------------------
#Sources
version_number_SOURCES = \
	version_number.f90


clean-local:
	rm -f $(top_builddir)/third_party_open/version_number/bin/linux/version_number
all-local: version_number
	${INSTALL} -d $(top_builddir)/third_party_open/version_number/bin/linux
	${INSTALL}    $(top_builddir)/third_party_open/version_number/packages/version_number/src/version_number $(top_builddir)/third_party_open/version_number/bin/linux

