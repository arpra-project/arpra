/*
 * arpra-test.h -- Header file for common testing routines.
 *
 * Copyright 2017-2018 James Paul Turner.
 *
 * This file is part of the Arpra library.
 *
 * The Arpra library is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published
 * by the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * The Arpra library is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
 * or FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public
 * License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with the Arpra library. If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef ARPRA_TEST_H
#define ARPRA_TEST_H

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif // HAVE_CONFIG_H

#include <stdlib.h>
#include <stdio.h>
#include <stdarg.h>
#include <string.h>
#include <time.h>

#include <arpra.h>
#ifdef WITH_MPFI
#include <arpra_to_mpfi.h>
#endif // WITH_MPFI

// RNG output mode enumeration.
enum test_rand_mode
{
    TEST_RAND_MIXED = 0,  // (-oo <  z  < +oo)
    TEST_RAND_SMALL_POS,  // (+0 <=  z  <  +1)
    TEST_RAND_SMALL_NEG,  // (-1  <  z  <= -0)
    TEST_RAND_LARGE_POS,  // (+1 <=  z  < +oo)
    TEST_RAND_LARGE_NEG,  // (-oo <  z  <= -1)
    TEST_RAND_SMALL,      // (+0  < |z| <  +1)
    TEST_RAND_LARGE,      // (+1 <= |z| < +oo)
    TEST_RAND_POS,        // (+0 <=  z  < +oo)
    TEST_RAND_NEG,        // (-oo <  z  <= -0)
};

#ifdef __cplusplus
extern "C" {
#endif

// Global test variables.
extern arpra_range x_A, y_A, z_A;
extern int test_fixture_ready;
#ifdef WITH_MPFI
extern mpfi_t x_I, y_I, z_I, z_AI;
extern arpra_mpfr rdiam_I, rdiam_AI, rdiam_diff;
#endif // WITH_MPFI

// Global RNG variables.
extern gmp_randstate_t test_randstate;
extern int test_rand_ready;

// Global logfile variables.
extern FILE *test_log;
extern int test_log_ready;

// Test fixture functions.
void test_fixture_init (arpra_precision prec, arpra_precision prec_internal);
void test_fixture_clear ();

// RNG functions.
void test_rand_init ();
void test_rand_clear ();
void test_rand_mpfr (arpra_mpfr *z, enum test_rand_mode mode);
void test_rand_arpra (arpra_range *z,
                      enum test_rand_mode mode_centre,
                      enum test_rand_mode mode_deviations);
#ifdef WITH_MPFI
void test_rand_mpfi (mpfi_ptr z,
                     enum test_rand_mode mode_low,
                     enum test_rand_mode mode_high);
#endif // WITH_MPFI

// Logfile functions.
void test_log_init (const char *test_name);
void test_log_clear ();
void test_log_printf (const char *format, ...);
void test_log_mpfr (const arpra_mpfr *x, const char *var_name);
#ifdef WITH_MPFI
void test_log_mpfi (mpfi_srcptr x, const char *var_name);
#endif // WITH_MPFI

// Symbol adjustments.
void test_share_all_syms (arpra_range *x, arpra_range *y);
void test_share_rand_syms (arpra_range *x, arpra_range *y);

// Test functions.
int test_compare_arpra (const arpra_range *x, const arpra_range *y);
#ifdef WITH_MPFI
void test_univariate (
    void (*f_arpra) (arpra_range *z, const arpra_range *x),
    int  (*f_mpfi) (mpfi_ptr z, mpfi_srcptr x));
void test_bivariate (
    void (*f_arpra) (arpra_range *z, const arpra_range *x, const arpra_range *y),
    int  (*f_mpfi) (mpfi_ptr z, mpfi_srcptr x, mpfi_srcptr y));
#endif // WITH_MPFI

#ifdef __cplusplus
}
#endif

#endif // ARPRA_TEST_H
