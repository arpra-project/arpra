/*
 * compare_bounds.c -- Test an MPFA result against upper and lower bounds.
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

mpfa_int_t test_compare_bounds (mpfa_srcptr x, mpfr_srcptr lo, mpfr_srcptr hi)
{
    mpfa_int_t inside;
    mpfa_prec_t prec;
    mpfr_t x_lo, x_hi;

    // Get lower and upper bounds of x.
    prec = mpfa_get_prec(x);
    mpfr_init2(x_lo, prec);
    mpfr_init2(x_hi, prec);
    mpfa_get_bounds(x_lo, x_hi, x);

    // Return nonzero if x does not include [lo, hi].
    inside = mpfr_greater_p(lo, x_lo) || mpfr_less_p(hi, x_hi);

    // Clear temp vars.
    mpfr_clear(x_lo);
    mpfr_clear(x_hi);
    return inside;
}
