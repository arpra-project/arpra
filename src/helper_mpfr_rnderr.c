/*
 * helper_mpfr_rnderr.c -- Compute MPFR functions with rounding error.
 *
 * Copyright 2017-2020 James Paul Turner.
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
 * This function assumes that the MPFR function which computed y returned
 * a nonzero ternary value, and thus y is inexact.
 *
 * Note: exponents of MPFR numbers are one greater than eponents of equivalent
 * IEEE-754 floating-point numbers, since MPFR significands are in [0.5, 1.0).
 */

void arpra_helper_mpfr_rnderr (mpfr_ptr err, mpfr_rnd_t rnd, mpfr_srcptr y)
{
    mpfr_t temp;
    mpfr_exp_t e;
    mpfr_prec_t p;
    arpra_prec internal_prec;

    // Initialise vars.
    internal_prec = arpra_get_internal_precision();
    mpfr_init2(temp, internal_prec);

    // Was y flushed to zero?
    if (mpfr_zero_p(y)) {
        // Rounding error is nextabove(0).
        mpfr_set_zero(temp, 1);
        mpfr_nextabove(temp);
    }
    else {
        e = mpfr_get_exp(y);
        p = mpfr_get_prec(y);

        // Nearest or directed rounding?
        if ((rnd == MPFR_RNDN) || (rnd == MPFR_RNDNA)) {
            // Rounding error is 0.5 ULP(y) = 2^(e-p-1).
            mpfr_set_si_2exp(temp, 1, (e - p - 1), MPFR_RNDU);
        }
        else {
            // Rounding error is ULP(y) = 2^(e-p).
            mpfr_set_si_2exp(temp, 1, (e - p), MPFR_RNDU);
        }
    }

    // Add rounding error to total.
    mpfr_add(err, err, temp, MPFR_RNDU);

    // Clear vars.
    mpfr_clear(temp);
}
