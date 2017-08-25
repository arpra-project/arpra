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

#include "mpfa.h"

#ifdef __cplusplus
extern "C" {
#endif

// Initialise and clear.
void mpfa_test_begin ();
void mpfa_test_end ();

// Compare known affine form result.
int mpfa_test_mpfa_1 (void (*mpfa_1) (mpfa_ptr z, mpfa_srcptr x),
                      mpfa_ptr expect, mpfa_srcptr x);
int mpfa_test_mpfa_2 (void (*mpfa_2) (mpfa_ptr z, mpfa_srcptr x, mpfa_srcptr y),
                      mpfa_ptr expect, mpfa_srcptr x, mpfa_srcptr y);
int mpfa_test_cmp_mpfa (mpfa_srcptr op1, mpfa_srcptr op2);

#ifdef WITH_MPFI

#include <mpfi.h>

// Compare MPFI interval result.
int mpfa_test_mpfi_1 (void (*mpfa_1) (mpfa_ptr z, mpfa_srcptr x),
                      void (*mpfi_1) (mpfi_ptr z, mpfi_srcptr x),
                      mpfa_srcptr x_a);
int mpfa_test_mpfi_2 (void (*mpfa_2) (mpfa_ptr z, mpfa_srcptr x, mpfa_srcptr y),
                      void (*mpfi_2) (mpfi_ptr z, mpfi_srcptr x, mpfi_srcptr y),
                      mpfa_srcptr x_a, mpfa_srcptr y_a);
int mpfa_test_cmp_mpfi (mpfa_srcptr op1, mpfi_srcptr op2);

#endif // WITH_MPFI

#ifdef __cplusplus
}
#endif

#endif // MPFA_TEST_H
