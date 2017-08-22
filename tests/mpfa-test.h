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
#ifdef WITH_MPFI
#include <mpfi.h>
#endif // WITH_MPFI

#ifdef __cplusplus
extern "C" {
#endif

// Initialise and clear.
void mpfa_test_begin ();
void mpfa_test_end ();

// Testing functions.
int mpfa_test_mpfa (void (*mpfa_fun) (mpfa_ptr z, mpfa_srcptr x, mpfa_srcptr y),
                    mpfa_ptr expected, mpfa_srcptr x, mpfa_srcptr y);
#ifdef WITH_MPFI
int mpfa_test_mpfi (void (*mpfa_fun) (mpfa_ptr z, mpfa_srcptr x, mpfa_srcptr y),
                    void (*mpfi_fun) (mpfi_ptr z, mpfi_srcptr x, mpfi_srcptr y),
                    mpfa_srcptr x, mpfa_srcptr y);
#endif // WITH_MPFI

#ifdef __cplusplus
}
#endif

#endif // MPFA_TEST_H
