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

# Programs
AC_PROG_CC
AC_PROG_CPP
AM_PROG_AR

LT_PREREQ([2.2])
LT_INIT([])

# GNU GMP
AC_CHECK_HEADER([gmp.h], [],
  [AC_MSG_ERROR([GMP header (gmp.h) is missing or unusable - see README])])
AC_CHECK_LIB([gmp], [__gmpz_init], [],
  [AC_MSG_ERROR([GMP library (libgmp.so) is missing or unusable - see README])])

# GNU MPFR
AC_CHECK_HEADER([mpfr.h], [],
  [AC_MSG_ERROR([MPFR header (mpfr.h) is missing or unusable - see README])])
AC_CHECK_LIB([mpfr], [mpfr_init], [],
  [AC_MSG_ERROR([MPFR library (libmpfr.so) is missing or unusable - see README])])

# Optionally enable MPFI tests (requires MPFI)
AC_ARG_ENABLE([mpfi-tests],
  [AS_HELP_STRING([--enable-mpfi-tests], [enable MPFI tests [default=no]])],
  [case "${enableval}" in
    yes) mpfi-tests=true ;;
    no)  mpfi-tests=false ;;
    *) AC_MSG_ERROR([bad value ${enableval} for --enable-mpfi-tests]) ;;
  esac
  # MPFI (required for MPFI tests)
  AC_CHECK_HEADER([mpfi.h], [],
    [AC_MSG_ERROR([MPFI header (mpfi.h) is missing or unusable - see README])])
  AC_CHECK_LIB([mpfi], [mpfi_init], [],
    [AC_MSG_ERROR([MPFI library (libmpfi.so) is missing or unusable - see README])])],
  [mpfi-tests=false])
AM_CONDITIONAL([MPFI_TESTS], [test x$mpfi-tests = xtrue])

# Header files
AC_CHECK_HEADERS([stdlib.h])

# Library functions
AC_CHECK_FUNCS([sqrt])
AC_FUNC_MALLOC
AC_FUNC_REALLOC

# Output
AC_CONFIG_FILES([Makefile])
AC_OUTPUT
