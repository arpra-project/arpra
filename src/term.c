/*
 * term.c -- Compute a single deviation term of an affine form.
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

#include "mpfa.h"
#include <assert.h>

mpfa_int_t mpfa_term (mpfr_ptr z, mpfr_srcptr x, mpfr_srcptr y, mpfr_srcptr alpha, mpfr_srcptr beta, mpfr_srcptr gamma) {
    mpfa_int_t inexact;
    mpfr_t alpha_x, beta_y;

    mpfr_init2(alpha_x, (mpfr_get_prec(alpha) + mpfr_get_prec(x)));
    assert(mpfr_mul(alpha_x, alpha, x, MPFR_RNDN) == 0);

    mpfr_init2(beta_y, (mpfr_get_prec(beta) + mpfr_get_prec(y)));
    assert(mpfr_mul(beta_y, beta, y, MPFR_RNDN) == 0);

    if (gamma == NULL) {
        inexact = mpfr_add(z, alpha_x, beta_y, MPFR_RNDN);
    }
    else {
        mpfr_ptr sumArray[3];
        sumArray[0] = alpha_x;
        sumArray[1] = beta_y;
        sumArray[2] = (mpfr_ptr) gamma;
        inexact = mpfr_sum(z, sumArray, 3, MPFR_RNDN);
    }

    mpfr_clear(alpha_x);
    mpfr_clear(beta_y);
    return inexact;
}
