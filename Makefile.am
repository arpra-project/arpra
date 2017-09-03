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

# Public headers
include_HEADERS = include/mpfa.h
if WITH_MPFI
include_HEADERS += include/mpfa2mpfi.h
endif WITH_MPFI

# MPFA library
lib_LTLIBRARIES = lib/libmpfa.la

lib_libmpfa_la_CPPFLAGS = -I$(top_srcdir)/include
lib_libmpfa_la_SOURCES = \
	src/mpfa-impl.h	src/add.c src/affine_1.c src/affine_2.c		\
	src/clear.c src/condense_last_n.c src/condense_small.c		\
	src/default_prec.c src/div.c src/error.c src/exp.c src/init.c	\
	src/init2.c src/internal_prec.c src/inv.c src/log.c src/mul.c	\
	src/neg.c src/next_sym.c src/prec.c src/set_d.c			\
	src/set_d_rad.c	src/set_str.c src/set_str_rad.c src/set_mpfr.c	\
	src/set_mpfr_rad.c src/set.c src/sqrt.c	src/sub.c src/sum.c	\
	src/term.c
if WITH_MPFI
lib_libmpfa_la_SOURCES += src/get_mpfi.c src/set_mpfi.c
endif WITH_MPFI

# Testsuite helper library
check_LTLIBRARIES = tests/libmpfa-test.la

tests_libmpfa_test_la_CPPFLAGS = -I$(top_srcdir)/include
tests_libmpfa_test_la_LIBADD = lib/libmpfa.la
tests_libmpfa_test_la_SOURCES = \
	tests/mpfa-test.h tests/rand.c tests/rand_ui.c			\
	tests/rand_mpfr.c tests/rand_mpfa.c tests/cmp_mpfa.c		\
	tests/share_syms.c
if WITH_MPFI
tests_libmpfa_test_la_SOURCES += tests/rand_mpfi.c tests/cmp_mpfi.c
endif WITH_MPFI

# Testsuite programs
check_PROGRAMS = tests/test1 tests/t_add tests/morris_lecar

tests_test1_CPPFLAGS = -I$(top_srcdir)/include
tests_test1_LDADD = tests/libmpfa-test.la
tests_test1_SOURCES = tests/test1.c

tests_t_add_CPPFLAGS = -I$(top_srcdir)/include
tests_t_add_LDADD = tests/libmpfa-test.la
tests_t_add_SOURCES = tests/t_add.c

tests_morris_lecar_CPPFLAGS = -I$(top_srcdir)/include
tests_morris_lecar_LDADD = tests/libmpfa-test.la
tests_morris_lecar_SOURCES = tests/morris_lecar.c

TESTS = $(check_PROGRAMS)

# Clean miscellaneous files
mostlyclean-local:
	-rm -f *.dat
