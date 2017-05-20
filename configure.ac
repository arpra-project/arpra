# Copyright 2016-2017 James Paul Turner.
#
# This file is part of the MPFA library.
# Process this file with autoconf to produce a configure script.
#
# The MPFA library is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published
# by the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# The MPFA library is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
# or FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public
# License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with the MPFA library. If not, see <http://www.gnu.org/licenses/>.

AC_PREREQ([2.65])
AC_INIT([MPFA], [1.0.0-dev], [https://github.com/jamesturner246/mpfa/issues])

AC_CONFIG_AUX_DIR([build-aux])
AC_CONFIG_SRCDIR([include/mpfa.h])
AC_CONFIG_HEADER([config.h])
AC_CONFIG_MACRO_DIR([m4])

AM_INIT_AUTOMAKE([1.13 -Wall -Werror subdir-objects dist-bzip2])

# Checks for programs.
AC_PROG_CC
AC_PROG_CPP
AM_PROG_AR

LT_PREREQ([2.2])
LT_INIT([])

# Checks for GNU GMP header and library.
AC_CHECK_HEADER([gmp.h], [],
  [AC_MSG_ERROR([Cannot include GMP header - see https://gmplib.org/])])
AC_CHECK_LIB([gmp], [__gmpz_init], [],
  [AC_MSG_ERROR([Cannot link GMP library - see https://gmplib.org/])])

# Checks for GNU MPFR header and library.
AC_CHECK_HEADER([mpfr.h], [],
  [AC_MSG_ERROR([Cannot include MPFR header - see http://www.mpfr.org/])])
AC_CHECK_LIB([mpfr], [mpfr_init], [],
  [AC_MSG_ERROR([Cannot link MPFR library - see http://www.mpfr.org/])])

# Checks for header files.
AC_CHECK_HEADERS([stdlib.h])

# Checks for library functions.
AC_FUNC_MALLOC
AC_FUNC_REALLOC

AC_CONFIG_FILES([Makefile])
AC_OUTPUT
