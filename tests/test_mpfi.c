/*
 * test_mpfi.c -- Test an affine form against an MPFI interval.
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

int mpfa_test_mpfi (mpfa_srcptr x, mpfi_srcptr y)
{
    mpfi_t x_i;
    mpfa_prec_t prec_x;

    // Convert x to an MPFI interval.
    prec_x = mpfa_get_prec(x);
    mpfi_init2(x_i, prec_x);
    mpfa_get_mpfi(x_i, x_a);

    // Return 1 if y is not a subinterval of x.
    if (mpfr_greater_p(&(x_i->left), &(y_i->left))) {
        mpfi_clear(x_i);
        return 1;
    }

    if (mpfr_less_p(&(x_i->right), &(y_i->right))) {
        mpfi_clear(x_i);
        return 1;
    }

    // Else return 0.
    mpfi_clear(x_i);
    return 0;
}

#endif // WITH_MPFI
