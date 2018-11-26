# Copyright 2016-2018 James Paul Turner.
#
# This file is part of the Arpra library.
# Process this file with autoconf to produce a configure script.
#
# The Arpra library is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published
# by the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# The Arpra library is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
# or FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public
# License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with the Arpra library. If not, see <http://www.gnu.org/licenses/>.

AC_PREREQ([2.65])
AC_INIT([Arpra], [0.1.0], [https://github.com/arpra-project/arpra/issues])

AC_CONFIG_AUX_DIR([build-aux])
AC_CONFIG_SRCDIR([include/arpra.h])
AC_CONFIG_HEADER([config.h])
AC_CONFIG_MACRO_DIR([m4])

AM_INIT_AUTOMAKE([1.13 -Wall -Werror subdir-objects dist-bzip2])

# Programs
AC_PROG_CC
AC_PROG_CPP
AM_PROG_AR

LT_PREREQ([2.2])
LT_INIT([])

# Standard libs
AC_CHECK_LIB([m], [ceil], [],
  [AC_MSG_ERROR([libm is missing or unusable])])
AC_CHECK_LIB([rt], [clock_gettime], [],
  [AC_MSG_ERROR([librt is missing or unusable])])

# GNU GMP
AC_CHECK_HEADER([gmp.h], [],
  [AC_MSG_ERROR([GMP header is missing or unusable - see README])])
AC_CHECK_LIB([gmp], [__gmpz_init], [],
  [AC_MSG_ERROR([GMP library is missing or unusable - see README])])

# GNU MPFR
AC_CHECK_HEADER([mpfr.h], [],
  [AC_MSG_ERROR([MPFR header is missing or unusable - see README])])
AC_CHECK_LIB([mpfr], [mpfr_init], [],
  [AC_MSG_ERROR([MPFR library is missing or unusable - see README])])

# MPFI
AC_CHECK_HEADER([mpfi.h], [],
  [AC_MSG_ERROR([MPFI header is missing or unusable - see README])])
AC_CHECK_LIB([mpfi], [mpfi_init], [],
  [AC_MSG_ERROR([MPFI library is missing or unusable - see README])])

# Header files
AC_CHECK_HEADERS([stdlib.h])

# Library functions
AC_CHECK_FUNCS([strtoul clock_gettime sqrt])
AC_FUNC_MALLOC
AC_FUNC_REALLOC

# Output
AC_CONFIG_FILES([Makefile])
AC_OUTPUT
