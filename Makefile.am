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
include_HEADERS = include/mpfa.h

lib_LTLIBRARIES = lib/libmpfa.la

lib_libmpfa_la_CPPFLAGS = -I$(top_srcdir)/include
lib_libmpfa_la_SOURCES = src/add.c \
                         src/affine_1.c \
                         src/affine_2.c \
                         src/clear.c \
                         src/condense_last_n.c \
                         src/condense_small.c \
                         src/default_prec.c \
                         src/div.c \
                         src/error.c \
                         src/exp.c \
                         src/get_prec.c \
                         src/init.c \
                         src/init2.c \
                         src/internal_prec.c \
                         src/inv.c \
                         src/log.c \
                         src/mul.c \
                         src/neg.c \
                         src/next_sym.c \
                         src/set_d.c \
                         src/set_d_rad.c \
                         src/set_mpfr.c \
                         src/set_mpfr_rad.c \
                         src/set_prec.c \
                         src/set.c \
                         src/sqrt.c \
                         src/sub.c \
                         src/sum.c \
                         src/term.c

# MPFA test suite
TESTS = tests/test1 tests/morris_lecar

check_PROGRAMS = tests/test1 tests/morris_lecar

tests_test1_CPPFLAGS = -I$(top_srcdir)/include
tests_test1_LDADD = lib/libmpfa.la
tests_test1_SOURCES = tests/test1.c

tests_morris_lecar_CPPFLAGS = -I$(top_srcdir)/include
tests_morris_lecar_LDADD = lib/libmpfa.la
tests_morris_lecar_SOURCES = tests/morris_lecar.c
