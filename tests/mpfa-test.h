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

// Global RNG variables.
extern gmp_randstate_t test_randstate;
extern int test_rand_ready;

// Global log variables.
extern FILE *test_log;
extern int test_log_ready;

// RNG initialise and clear.
void test_rand_init ();
void test_rand_clear ();

// Logfile initialise and clear.
void test_log_init (const char *test_name);
void test_log_clear ();

// Generate random arguments.
mpfa_uint_t test_rand_ui (mpfa_uint_t n);
void test_rand_mpfr (mpfr_ptr z, enum test_rand_mode mode);
void test_rand_mpfa (mpfa_ptr z, enum test_rand_mode mode);
#ifdef WITH_MPFI
void test_rand_mpfi (mpfi_ptr z, enum test_rand_mode mode);
#endif // WITH_MPFI

// Adjust affine form symbols.
void test_share_syms (mpfa_ptr x, mpfa_ptr y, const mpfa_uint_t share_chance);

// Compare functions.
mpfa_int_t test_compare_bounds (mpfa_srcptr x, mpfr_srcptr lo, mpfr_srcptr hi);
mpfa_int_t test_compare_mpfa (mpfa_srcptr x, mpfa_srcptr y);
#ifdef WITH_MPFI
mpfa_int_t test_compare_mpfi (mpfa_srcptr x, mpfi_srcptr y);
#endif // WITH_MPFI

#ifdef __cplusplus
}
#endif

#endif // MPFA_TEST_H
