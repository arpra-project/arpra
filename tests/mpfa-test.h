/*
 * mpfa-test.h -- Header file for common testing routines.
 *
 * Copyright 2017 James Paul Turner.
 *
 * This file is part of the MPFA library.
 *
 * The MPFA library is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published
 * by the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * The MPFA library is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
 * or FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public
 * License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with the MPFA library. If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef MPFA_TEST_H
#define MPFA_TEST_H

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif // HAVE_CONFIG_H

#include <stdlib.h>
#include <stdio.h>
#include <stdarg.h>
#include <string.h>
#include <time.h>

#include <mpfa.h>
#ifdef WITH_MPFI
#include <mpfa2mpfi.h>
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
extern mpfa_t x_A, y_A, z_A;
extern int test_fixture_ready;
#ifdef WITH_MPFI
extern mpfi_t x_I, y_I, z_I, z_AI;
extern mpfr_t rdiam_I, rdiam_AI, rdiam_diff;
#endif // WITH_MPFI

// Global RNG variables.
extern gmp_randstate_t test_randstate;
extern int test_rand_ready;

// Global logfile variables.
extern FILE *test_log;
extern int test_log_ready;

// Test fixture functions.
void test_fixture_init (mpfa_prec_t prec, mpfa_prec_t prec_internal);
void test_fixture_clear ();

// RNG functions.
void test_rand_init ();
void test_rand_clear ();
void test_rand_mpfr (mpfr_ptr z, enum test_rand_mode mode);
void test_rand_mpfa (mpfa_ptr z, enum test_rand_mode mode);
#ifdef WITH_MPFI
void test_rand_mpfi (mpfi_ptr z, enum test_rand_mode mode);
#endif // WITH_MPFI

// Logfile functions.
void test_log_init (const char *test_name);
void test_log_clear ();
void test_log_printf (const char *format, ...);
void test_log_mpfr (mpfr_srcptr x, const char *var_name);
#ifdef WITH_MPFI
void test_log_mpfi (mpfi_srcptr x, const char *var_name);
#endif // WITH_MPFI

// Symbol adjustments.
void test_share_all_syms (mpfa_ptr x, mpfa_ptr y);
void test_share_rand_syms (mpfa_ptr x, mpfa_ptr y);

// Compare functions.
int test_compare_mpfa (mpfa_srcptr x, mpfa_srcptr y);
#ifdef WITH_MPFI
int test_univariate_mpfi (
    void (*f_MPFA) (mpfa_ptr z, mpfa_srcptr x),
    int  (*f_MPFI) (mpfi_ptr z, mpfi_srcptr x));
int test_bivariate_mpfi (
    void (*f_MPFA) (mpfa_ptr z, mpfa_srcptr x, mpfa_srcptr y),
    int  (*f_MPFI) (mpfi_ptr z, mpfi_srcptr x, mpfi_srcptr y));
#endif // WITH_MPFI

#ifdef __cplusplus
}
#endif

#endif // MPFA_TEST_H
