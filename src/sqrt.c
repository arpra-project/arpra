/*
 * sqrt.c -- Compute the square root of an affine form.
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
 * This affine square root function uses a Chebyshev linear approximation.
 */

void mpfa_sqrt (mpfa_ptr z, mpfa_srcptr x)
{
    mpfr_t temp, xa, xb, da, db, du, alpha, gamma, delta;
    mpfa_prec_t prec_internal;

    // Init temp vars with internal precision.
    prec_internal = mpfa_get_internal_prec();
    mpfr_init2(temp, prec_internal);
    mpfr_init2(xa, prec_internal);
    mpfr_init2(xb, prec_internal);
    mpfr_init2(da, prec_internal);
    mpfr_init2(db, prec_internal);
    mpfr_init2(du, prec_internal);
    mpfr_init2(alpha, prec_internal);
    mpfr_init2(gamma, prec_internal);
    mpfr_init2(delta, prec_internal);

    // Handle trivial case, where x has zero radius.
    if (mpfr_zero_p(&(x->radius))) {
        if (mpfr_sqrt(temp, &(x->centre), MPFR_RNDN)) {
            mpfa_error(delta, temp);
        }
        else {
            mpfr_set_si(delta, 0, MPFR_RNDN);
        }

        mpfa_set_mpfr_rad(z, temp, delta);
    }
    else {
        mpfr_sub(xa, &(x->centre), &(x->radius), MPFR_RNDD);
        mpfr_add(xb, &(x->centre), &(x->radius), MPFR_RNDU);

        // Set NaN if x contains some p <= 0.
        if (mpfr_sgn(xa) <= 0) {
            mpfa_clear_terms(z);

            // TODO: find a better representation for Inf
            mpfr_set_nan(&(z->centre));
        }
        else {
            // compute alpha
            mpfr_sqrt(alpha, xa, MPFR_RNDN);
            mpfr_sqrt(temp, xb, MPFR_RNDN);
            mpfr_add(alpha, alpha, temp, MPFR_RNDN);
            mpfr_si_div(alpha, 1, alpha, MPFR_RNDN);

            // compute difference (sqrt(a) - alpha a)
            mpfr_mul(da, alpha, xa, MPFR_RNDU);
            mpfr_sqrt(temp, xa, MPFR_RNDD);
            mpfr_sub(da, temp, da, MPFR_RNDD);

            // compute difference (sqrt(b) - alpha b)
            mpfr_mul(db, alpha, xb, MPFR_RNDU);
            mpfr_sqrt(temp, xb, MPFR_RNDD);
            mpfr_sub(db, temp, db, MPFR_RNDD);

            mpfr_min(da, da, db, MPFR_RNDN);

            // compute difference (sqrt(u) - alpha u)
            mpfr_si_div(du, 1, alpha, MPFR_RNDU);
            mpfr_div_si(du, du, 4, MPFR_RNDU);

            // compute gamma
            mpfr_add(gamma, da, du, MPFR_RNDN);
            mpfr_div_si(gamma, gamma, 2, MPFR_RNDN);

            // compute delta
            mpfr_sub(delta, du, gamma, MPFR_RNDU);
            mpfr_sub(temp, gamma, da, MPFR_RNDU);
            mpfr_max(delta, delta, temp, MPFR_RNDU);

            // compute affine approximation
            mpfa_affine_1(z, x, alpha, gamma, delta);
        }
    }

    // Clear temp vars.
    mpfr_clear(temp);
    mpfr_clear(xa);
    mpfr_clear(xb);
    mpfr_clear(da);
    mpfr_clear(db);
    mpfr_clear(du);
    mpfr_clear(alpha);
    mpfr_clear(gamma);
    mpfr_clear(delta);
}
