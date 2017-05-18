# Copyright 2016-2017 James Paul Turner.
#
# This file is part of the MPFA library.
# Process this file with automake to produce a Makefile.in file.
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

ACLOCAL_AMFLAGS = -I m4

# MPFA library
include_HEADERS = src/mpfa.h

lib_LTLIBRARIES = libmpfa.la

libmpfa_la_SOURCES = src/add.c \
		     src/affine_1.c \
		     src/affine_2.c \
		     src/clear.c \
		     src/div.c \
		     src/exp.c \
		     src/get_default_prec.c \
		     src/get_prec.c \
		     src/init.c \
		     src/init2.c \
		     src/inv.c \
		     src/mul.c \
		     src/neg.c \
		     src/next_sym.c \
		     src/set_d.c \
		     src/set_default_prec.c \
		     src/set_mpfr.c \
		     src/set_prec.c \
		     src/set.c \
		     src/sub.c \
		     src/sum.c \
		     src/term.c

# MPFA test suite
TESTS = tests/mpfa_test1

check_PROGRAMS = tests/mpfa_test1

tests_mpfa_test1_SOURCES = tests/mpfa_test1.c
tests_mpfa_test1_CPPFLAGS = -I$(srcdir)/src
tests_mpfa_test1_LDADD = libmpfa.la

# MPFA documentation
