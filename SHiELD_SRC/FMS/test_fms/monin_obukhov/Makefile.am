#***********************************************************************
#*                   GNU Lesser General Public License
#*
#* This file is part of the GFDL Flexible Modeling System (FMS).
#*
#* FMS is free software: you can redistribute it and/or modify it under
#* the terms of the GNU Lesser General Public License as published by
#* the Free Software Foundation, either version 3 of the License, or (at
#* your option) any later version.
#*
#* FMS is distributed in the hope that it will be useful, but WITHOUT
#* ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
#* FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
#* for more details.
#*
#* You should have received a copy of the GNU Lesser General Public
#* License along with FMS.  If not, see <http://www.gnu.org/licenses/>.
#***********************************************************************

# This is an automake file for the test_fms/monin_obukhov directory of the FMS
# package.

# uramirez, Ed Hartnett

# Find the needed mod files.
AM_CPPFLAGS = -I$(MODDIR)

# Link to the FMS library.
LDADD = $(top_builddir)/libFMS/libFMS.la

# Build this test program.
check_PROGRAMS = test_monin_obukhov

# This is the source code for the test.
test_monin_obukhov_SOURCES = test_monin_obukhov.F90

TEST_EXTENSIONS = .sh
SH_LOG_DRIVER = env AM_TAP_AWK='$(AWK)' $(SHELL) \
                  $(top_srcdir)/tap-driver.sh

# Run the test program.
TESTS = test_monin_obukhov2.sh

# These files will also be included in the distribution.
EXTRA_DIST = test_monin_obukhov2.sh

# Clean up
CLEANFILES = input.nml *.out *.dpi *.spi *.dyn *.spl
