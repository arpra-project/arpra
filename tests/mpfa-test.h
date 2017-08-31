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

#include <mpfa.h>
#ifdef WITH_MPFI
#include <mpfa2mpfi.h>
#endif // WITH_MPFI

// Modes of RNG output.
enum mpfa_test_rand_mode {
    MPFA_TEST_RAND_MIXED,
    MPFA_TEST_RAND_SMALL,
    MPFA_TEST_RAND_LARGE,
    MPFA_TEST_RAND_POS,
    MPFA_TEST_RAND_NEG
};

#ifdef __cplusplus
extern "C" {
#endif

// Initialise and clear RNG.
void mpfa_test_rand_init ();
void mpfa_test_rand_clear ();

// Generate random arguments.
mpfa_uint_t mpfa_test_rand_ui (mpfa_uint_t n_bits);
void mpfa_test_rand_mpfr (mpfr_ptr z, enum mpfa_test_rand_mode mode);
void mpfa_test_rand_mpfa (mpfa_ptr z, enum mpfa_test_rand_mode mode);
#ifdef WITH_MPFI
void mpfa_test_rand_mpfi (mpfi_ptr z, enum mpfa_test_rand_mode mode);
#endif // WITH_MPFI

// Compare function results.
mpfa_int_t mpfa_test_cmp_mpfa (mpfa_srcptr x, mpfa_srcptr y);
#ifdef WITH_MPFI
mpfa_int_t mpfa_test_cmp_mpfi (mpfa_srcptr x, mpfi_srcptr y);
#endif // WITH_MPFI

#ifdef __cplusplus
}
#endif

#endif // MPFA_TEST_H
