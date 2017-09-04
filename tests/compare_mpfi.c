/*
 * compare_mpfi.c -- Test an MPFA result against an MPFI result.
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

#include "mpfa-test.h"

mpfa_int_t test_compare_mpfi (mpfa_srcptr x, mpfi_srcptr y)
{
    mpfa_int_t inside;
    mpfa_prec_t prec_x;
    mpfi_t x_i;

    // Convert x to an MPFI interval.
    prec_x = mpfa_get_prec(x);
    mpfi_init2(x_i, prec_x);
    mpfa_get_mpfi(x_i, x);

    // Return 1 if x is inside y, else 0.
    inside = mpfi_is_strictly_inside(x_i, y);

    // Clear temp vars.
    mpfi_clear(x_i);
    return inside;
}
