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
#ifdef WITH_MPFI
#include <mpfi.h>
#endif // WITH_MPFI
#include "mpfa.h"

#ifdef __cplusplus
extern "C" {
#endif

// RNG initialise and clear.
void mpfa_test_rand_init ();
void mpfa_test_rand_clear ();

// Test functions.
int mpfa_test_mpfa (mpfa_srcptr x, mpfa_srcptr y);

#ifdef WITH_MPFI

// MPFI test functions.
int mpfa_test_mpfi (mpfa_srcptr x, mpfi_srcptr y);

#endif // WITH_MPFI

#ifdef __cplusplus
}
#endif

#endif // MPFA_TEST_H
