/*
 * helper_error.c -- Adds rounding error from inexact MPFR operations.
 *
 * Copyright 2017-2018 James Paul Turner.
 *
 * This file is part of the Arpra library.
 *
 * The Arpra library is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published
 * by the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * The Arpra library is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
 * or FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public
 * License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with the Arpra library. If not, see <http://www.gnu.org/licenses/>.
 */

#include "arpra-impl.h"

/*
 * These functions assume that x is inexact - i.e. the last MPFR function
 * of x returned a nonzero integer.
 *
 * MPFR exponents are one less than exponents of equivalent IEEE-754
 * floating-point numbers, since MPFR significands are in [0.5, 1.0).
 */

/*
 * arpra_helper_error_ulp:
 * If x DID NOT underflow, then added error is ULP(x).
 * If x DID underflow, then added error is nextabove(0).
 */

void arpra_helper_error_ulp (arpra_mpfr *error, const arpra_mpfr *x)
{
    arpra_mpfr temp;
    arpra_prec prec_internal;
    mpfr_prec_t p;
    mpfr_exp_t e;

    // Initialise vars.
    prec_internal = arpra_get_internal_precision();
    mpfr_init2(&temp, prec_internal);

    if (mpfr_zero_p(x)) {
        // Added error is nextabove(0)
        mpfr_set_zero(&temp, 1);
        mpfr_nextabove(&temp);
    }
    else {
        // Added error is ULP(x) = 2^(e-p)
        p = mpfr_get_prec(x);
        e = mpfr_get_exp(x);
        mpfr_set_si_2exp(&temp, 1, (e - p), MPFR_RNDU);
    }

    // Add error to total.
    mpfr_add(error, error, &temp, MPFR_RNDU);

    // Clear vars.
    mpfr_clear(&temp);
}

/*
 * arpra_helper_error_half_ulp:
 * If x DID NOT underflow, then added error is 1/2 ULP(x).
 * If x DID underflow, then added error is nextabove(0).
 */

void arpra_helper_error_half_ulp (arpra_mpfr *error, const arpra_mpfr *x)
{
    arpra_mpfr temp;
    arpra_prec prec_internal;
    mpfr_prec_t p;
    mpfr_exp_t e;

    // Initialise vars.
    prec_internal = arpra_get_internal_precision();
    mpfr_init2(&temp, prec_internal);

    if (mpfr_zero_p(x)) {
        // Added error is nextabove(0)
        mpfr_set_zero(&temp, 1);
        mpfr_nextabove(&temp);
    }
    else {
        // Added error is 1/2 ULP(x) = 2^(e-p-1)
        p = mpfr_get_prec(x);
        e = mpfr_get_exp(x);
        mpfr_set_si_2exp(&temp, 1, (e - p - 1), MPFR_RNDU);
    }

    // Add error to total.
    mpfr_add(error, error, &temp, MPFR_RNDU);

    // Clear vars.
    mpfr_clear(&temp);
}
