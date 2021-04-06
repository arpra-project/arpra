# Copyright 2016-2020 James Paul Turner.
#
# This file is part of the Arpra library.
# Process this file with automake to produce a Makefile.in file.
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

ACLOCAL_AMFLAGS = -I m4
AM_CPPFLAGS = -I$(top_srcdir)/include

# Public headers
include_HEADERS = include/arpra.h include/arpra_ode.h

# Arpra library
lib_LTLIBRARIES = lib/libarpra.la
lib_libarpra_la_SOURCES = \
	src/arpra-impl.h src/add.c src/sub.c src/helper_affine_1.c	\
	src/helper_affine_2.c src/clear.c src/reduce_last_n.c		\
	src/reduce_small.c src/default_precision.c src/div.c src/exp.c	\
	src/predicates.c src/init.c src/internal_precision.c src/inv.c	\
	src/log.c src/mul.c src/neg.c src/helper_symbol.c		\
	src/precision.c src/set.c src/set_special.c src/sqrt.c		\
	src/sum.c src/get_bounds.c src/ode_stepper.c src/increase.c	\
	src/ode_euler.c src/ode_trapezoidal.c src/ode_bogsham32.c	\
	src/ode_dopri54.c src/ode_dopri87.c src/helper_mpfr_rnderr.c	\
	src/helper_buffer.c src/ext_mpfr.c src/get_mpfi.c		\
	src/helper_compute_range.c src/helper_check_result.c		\
	src/set_mpfi.c src/mpfr_fn.c src/helper_clear_terms.c		\
	src/helper_mix_trim.c src/range_method.c

# Testsuite helper library
check_LTLIBRARIES = tests/libarpra-test.la
tests_libarpra_test_la_LIBADD = lib/libarpra.la
tests_libarpra_test_la_SOURCES = \
	tests/arpra-test.h tests/share_syms.c tests/fixture.c		\
	tests/logfile.c tests/logfile_printf.c tests/logfile_mpfr.c	\
	tests/rand.c tests/rand_mpfr.c tests/rand_arpra.c		\
	tests/compare_arpra.c tests/univariate.c tests/bivariate.c	\
	tests/logfile_mpfi.c

# Testsuite test programs
check_PROGRAMS = \
	tests/t_add tests/t_sub tests/t_mul tests/t_div	tests/t_neg	\
	tests/t_inv tests/t_sqrt tests/t_exp tests/t_log
tests_t_add_LDADD = tests/libarpra-test.la
tests_t_add_SOURCES = tests/t_add.c
tests_t_sub_LDADD = tests/libarpra-test.la
tests_t_sub_SOURCES = tests/t_sub.c
tests_t_mul_LDADD = tests/libarpra-test.la
tests_t_mul_SOURCES = tests/t_mul.c
tests_t_div_LDADD = tests/libarpra-test.la
tests_t_div_SOURCES = tests/t_div.c
tests_t_neg_LDADD = tests/libarpra-test.la
tests_t_neg_SOURCES = tests/t_neg.c
tests_t_inv_LDADD = tests/libarpra-test.la
tests_t_inv_SOURCES = tests/t_inv.c
tests_t_sqrt_LDADD = tests/libarpra-test.la
tests_t_sqrt_SOURCES = tests/t_sqrt.c
tests_t_exp_LDADD = tests/libarpra-test.la
tests_t_exp_SOURCES = tests/t_exp.c
tests_t_log_LDADD = tests/libarpra-test.la
tests_t_log_SOURCES = tests/t_log.c
TESTS = $(check_PROGRAMS)

# Extra programs
EXTRA_PROGRAMS =

EXTRA_PROGRAMS += extra/henon_map
extra_henon_map_LDADD = lib/libarpra.la
extra_henon_map_SOURCES = extra/henon_map.c

EXTRA_PROGRAMS += extra/henon_map_mpfi
extra_henon_map_mpfi_SOURCES = extra/henon_map_mpfi.c

EXTRA_PROGRAMS += extra/morris_lecar
extra_morris_lecar_LDADD = lib/libarpra.la
extra_morris_lecar_SOURCES = extra/morris_lecar.c

EXTRA_PROGRAMS += extra/morris_lecar_mpfr
extra_morris_lecar_mpfr_SOURCES = extra/morris_lecar_mpfr.c

EXTRA_PROGRAMS += extra/morris_lecar_mpfi
extra_morris_lecar_mpfi_SOURCES = extra/morris_lecar_mpfi.c

EXTRA_PROGRAMS += extra/hodgkin_huxley
extra_hodgkin_huxley_LDADD = lib/libarpra.la
extra_hodgkin_huxley_SOURCES = extra/hodgkin_huxley.c

# Experiments
EXTRA_PROGRAMS += experiments/experiment_1
experiments_experiment_1_LDADD = lib/libarpra.la
experiments_experiment_1_SOURCES = experiments/experiment_1.c

EXTRA_PROGRAMS += experiments/experiment_2
experiments_experiment_2_SOURCES = experiments/experiment_2.c

EXTRA_PROGRAMS += experiments/experiment_3/morris_lecar
experiments_experiment_3_morris_lecar_LDADD = lib/libarpra.la
experiments_experiment_3_morris_lecar_SOURCES = experiments/experiment_3/morris_lecar.c

EXTRA_PROGRAMS += experiments/fig5/ml-arpra
experiments_fig5_ml_arpra_LDADD = lib/libarpra.la
experiments_fig5_ml_arpra_SOURCES = experiments/fig5/ml-arpra.c

EXTRA_PROGRAMS += experiments/fig5/ml-mpfr
experiments_fig5_ml_mpfr_SOURCES = experiments/fig5/ml-mpfr.c

EXTRA_PROGRAMS += experiments/fig5/ml-mpfi
experiments_fig5_ml_mpfi_SOURCES = experiments/fig5/ml-mpfi.c

# Documentation
info_TEXINFOS = doc/arpra.texi
doc_arpra_TEXINFOS = doc/fdl-1.3.texi

# Clean miscellaneous files
mostlyclean-local:
	-rm -f *.log
