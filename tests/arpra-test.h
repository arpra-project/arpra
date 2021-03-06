/*
 * arpra-test.h -- Header file for common testing routines.
 *
 * Copyright 2017-2020 James Paul Turner.
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

#include "../src/arpra-impl.h"

// RNG output mode enumeration.
typedef enum test_rand_mode_enum test_rand_mode;
enum test_rand_mode_enum
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
extern int test_fixture_ready;
extern arpra_range y_A, x1_A, x2_A;
extern mpfi_t y_I;
extern mpfr_t y_I_diam, y_A_diam, y_A_diam_rel;

// Global RNG variables.
extern int test_rand_ready;
extern gmp_randstate_t test_randstate;

// Global logfile variables.
extern int test_log_ready;
extern FILE *test_log;

// Test fixture functions.
void test_fixture_init (arpra_prec prec, arpra_prec prec_internal);
void test_fixture_clear ();

// RNG functions.
void test_rand_init ();
void test_rand_clear ();
void test_rand_mpfr (mpfr_ptr y, arpra_prec prec, test_rand_mode mode);
void test_rand_arpra (arpra_range *y, test_rand_mode mode_c, test_rand_mode mode_d);
void test_rand_uniform_mpfr (mpfr_ptr y, long int y_a, long int y_b);
void test_rand_uniform_arpra (arpra_range *y,
                              long int yc_a, long int yc_b,
                              long int yd_a, long int yd_b);

// Logfile functions.
void test_log_init (const char *test_name);
void test_log_clear ();
void test_log_printf (const char *format, ...);
void test_log_mpfr (mpfr_srcptr x1, const char *var_name);
void test_log_mpfi (mpfi_srcptr x1, const char *var_name);

// Symbol adjustments.
void test_share_all_syms (arpra_range *x1, arpra_range *x2);
void test_share_rand_syms (arpra_range *x1, arpra_range *x2);
void test_share_n_syms (arpra_range *x1, arpra_range *x2, arpra_uint n);

// Test functions.
int test_compare_arpra (const arpra_range *x1, const arpra_range *x2);
void test_univariate (
    void (*f_arpra) (arpra_range *y, const arpra_range *x1),
    int  (*f_mpfi) (mpfi_ptr y, mpfi_srcptr x1));
void test_bivariate (
    void (*f_arpra) (arpra_range *y, const arpra_range *x1, const arpra_range *x2),
    int  (*f_mpfi) (mpfi_ptr y, mpfi_srcptr x1, mpfi_srcptr x2));

#ifdef __cplusplus
}
#endif

#endif // ARPRA_TEST_H
