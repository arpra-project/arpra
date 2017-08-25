/*
 * test_mpfi.c -- Compare function results against known MPFI intervals.
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

#ifdef WITH_MPFI

#include "mpfa-test.h"

int mpfa_test_mpfi_1 (void (*mpfa_1) (mpfa_ptr z, mpfa_srcptr x),
                      void (*mpfi_1) (mpfi_ptr z, mpfi_srcptr x),
                      mpfa_srcptr x_a)
{
    int success;
    mpfa_prec_t prec;
    mpfa_t z_a;
    mpfi_t z_i, x_i;

    // Init test vars with precision of input affine form.
    prec = mpfa_get_prec(x_a);
    mpfa_init2(z_a, prec);
    mpfi_init2(z_i, prec);
    mpfi_init2(x_i, prec);

    //mpfa_get_mpfi(x_i, x_a);
    //mpfa_1(z_a, x_a);
    //mpfi_1(z_i, x_i);
    //success = mpfa_test_cmp_mpfi(z_a, z_i);

    // Clear test vars and return test result.
    mpfa_clear(z_a);
    mpfi_clear(z_i);
    mpfi_clear(x_i);
    return success;
}

int mpfa_test_mpfi_2 (void (*mpfa_2) (mpfa_ptr z, mpfa_srcptr x, mpfa_srcptr y),
                      void (*mpfi_2) (mpfi_ptr z, mpfi_srcptr x, mpfi_srcptr y),
                      mpfa_srcptr x_a, mpfa_srcptr y_a)
{
    return 0;
}

int mpfa_test_cmp_mpfi (mpfa_srcptr op1, mpfi_srcptr op2)
{
    return 0;
}

#endif // WITH_MPFI
