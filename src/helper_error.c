/*
 * helper_error.c -- Compute the error of an inexact MPFR operation.
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
 * This function assumes that rounding error has occured - i.e. an MPFR call
 * has returned nonzero. If x == 0, then nextabove(x) is used instead of
 * 1/2 ULP(x). This ensures that error is always recorded, even when x was
 * flushed to zero.
 */

void arpra_helper_error (arpra_mpfr *error, const arpra_mpfr *x)
{
    arpra_precision p;
    arpra_exponent e;

    if (mpfr_zero_p(x)) {
        // error = nextabove(x)
        mpfr_set_si(error, 0, MPFR_RNDN);
        mpfr_nextabove(error);
    }
    else {
        // error = 1/2 ULP(x) = 2^(e-p-1)
        p = mpfr_get_prec(x);
        e = mpfr_get_exp(x);
        mpfr_set_si_2exp(error, 1, (e - p - 1), MPFR_RNDU);
    }
}
