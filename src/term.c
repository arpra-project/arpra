/*
 * term.c -- Compute a single deviation term.
 *
 * Copyright 2017 James Paul Turner.
 *
 * This file is part of the ArPRA library.
 *
 * The ArPRA library is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published
 * by the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * The ArPRA library is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
 * or FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public
 * License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with the ArPRA library. If not, see <http://www.gnu.org/licenses/>.
 */

#include "arpra-impl.h"

arpra_int_t arpra_term (mpfr_ptr z, mpfr_srcptr x, mpfr_srcptr y, mpfr_srcptr alpha, mpfr_srcptr beta, mpfr_srcptr gamma)
{
    arpra_int_t inexact;
    mpfr_t alpha_x, beta_y;

    // alpha * x needs precision prec(alpha) + prec(x) to be exact.
    mpfr_init2(alpha_x, (mpfr_get_prec(alpha) + mpfr_get_prec(x)));
    mpfr_mul(alpha_x, alpha, x, MPFR_RNDN);

    // beta * y needs precision prec(beta) + prec(y) to be exact.
    mpfr_init2(beta_y, (mpfr_get_prec(beta) + mpfr_get_prec(y)));
    mpfr_mul(beta_y, beta, y, MPFR_RNDN);

    if (gamma == NULL) {
        // z = (alpha * x) + (beta * y)
        inexact = mpfr_add(z, alpha_x, beta_y, MPFR_RNDN);
    }
    else {
        // z = (alpha * x) + (beta * y) + gamma
        mpfr_ptr sumArray[3];
        sumArray[0] = alpha_x;
        sumArray[1] = beta_y;
        sumArray[2] = (mpfr_ptr) gamma;
        inexact = mpfr_sum(z, sumArray, 3, MPFR_RNDN);
    }

    // Clear temp vars.
    mpfr_clear(alpha_x);
    mpfr_clear(beta_y);
    return inexact;
}
