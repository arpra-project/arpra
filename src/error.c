/*
 * error.c -- Compute the max ULP error of an inexact MPFR operation.
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

#include "mpfa-impl.h"

/*
 * This function assumes that rounding error has occured - i.e. an MPFR call
 * has returned nonzero. Therefore, if x == 0, then nextabove(x) is used
 * instead of 1/2 ULP(x). This ensures that error is always recorded, even
 * when x was flushed to zero.
 */

void mpfa_error (mpfr_ptr error, mpfr_srcptr x) {
    mpfa_prec_t p;
    mpfa_exp_t e;

    if (mpfr_zero_p(x)) {
        // error = nextabove(x)
        mpfr_set_si(error, 0, MPFR_RNDN);
        mpfr_nextabove(error);
    }
    else {
        // error = 1/2 ULP(x) = 2^(e-p-1)
        p = mpfr_get_prec(x);
        e = mpfr_get_exp(x);
        mpfr_set_si(error, (e - p - 1), MPFR_RNDU);
        mpfr_exp2(error, error, MPFR_RNDU);
    }
}
