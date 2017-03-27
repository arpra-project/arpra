#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

# Initialisation.
AC_PREREQ([2.65])
AC_COPYRIGHT([
Copyright 2016-2017 James Paul Turner.

This file is part of the MPFA library.

The MPFA library is free software: you can redistribute it and/or modify
it under the terms of the GNU Lesser General Public License as published
by the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

The MPFA library is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
or FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public
License for more details.

You should have received a copy of the GNU Lesser General Public License
along with the MPFA library. If not, see <http://www.gnu.org/licenses/>.
])

AC_CONFIG_SRCDIR([next_sym.c])
AC_CONFIG_HEADERS([config.h])
AC_CONFIG_FILES([Makefile src/Makefile])
AC_CONFIG_MACRO_DIR([m4])
AC_CONFIG_AUX_DIR([build-aux])

AC_INIT([MPFA], [1.0.0-dev])
AM_INIT_AUTOMAKE([-Wall -Werror])

# Checks for programs.
AC_PROG_CC

# Checks for libraries.
AC_CHECK_LIB([gmp], [__gmpz_init])
AC_CHECK_LIB([mpfr], [mpfr_init])
#AC_CHECK_LIB([gmp], [__gmpz_init], ,
#  [AC_MSG_ERROR([GNU MP not found, see https://gmplib.org/])])

# Checks for header files.
AC_CHECK_HEADERS([stdlib.h stdio.h assert.h gmp.h mpfr.h])
#AC_CHECK_HEADER([stdio.h])
#AC_CHECK_HEADER([stdlib.h])
#AC_CHECK_HEADER([assert.h])
#AC_CHECK_HEADER([gmp.h])
#AC_CHECK_HEADER([mpfr.h])

# Checks for typedefs, structures, and compiler characteristics.

# Checks for library functions.
AC_FUNC_MALLOC
AC_FUNC_REALLOC

AC_OUTPUT
