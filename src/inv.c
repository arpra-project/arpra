/*
 * inv.c -- Compute the inverse of an Arpra range.
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
 * This affine inverse function uses a Chebyshev linear approximation.
 */

void arpra_inv (arpra_range *z, const arpra_range *x)
{
    arpra_mpfr temp, xa, xb, da, db, du, alpha, gamma, delta;
    arpra_precision prec_internal;
    arpra_int sign;

    // Initialise vars.
    prec_internal = arpra_get_internal_precision();
    mpfr_init2(&temp, prec_internal);
    mpfr_init2(&xa, prec_internal);
    mpfr_init2(&xb, prec_internal);
    mpfr_init2(&da, prec_internal);
    mpfr_init2(&db, prec_internal);
    mpfr_init2(&du, prec_internal);
    mpfr_init2(&alpha, prec_internal);
    mpfr_init2(&gamma, prec_internal);
    mpfr_init2(&delta, prec_internal);

    // Handle x with zero radius.
    if (mpfr_zero_p(&(x->radius))) {
        if (mpfr_ui_div(&temp, 1, &(x->centre), MPFR_RNDN)) {
            arpra_helper_error(&delta, &temp);
            arpra_set_mpfr_rad(z, &temp, &delta);
        }
        else {
            arpra_set_mpfr(z, &temp);
        }
    }
    else {
        // Handle domain violations.
        if (arpra_nan_p(x)) {
            arpra_set_nan(z);
        }
        else if (arpra_has_zero_p(x)) {
            arpra_set_inf(z);
        }

        // Domain is OK.
        else {
            arpra_get_bounds(&xa, &xb, x);

            // Handle negative x.
            sign = mpfr_sgn(&xa);
            if (sign < 0) {
                mpfr_set(&temp, &xa, MPFR_RNDN);
                mpfr_neg(&xa, &xb, MPFR_RNDN);
                mpfr_neg(&xb, &temp, MPFR_RNDN);
            }

            // compute alpha
            mpfr_si_div(&alpha, -1, &xb, MPFR_RNDN);
            mpfr_div(&alpha, &alpha, &xa, MPFR_RNDN);

            // compute difference (1/a - alpha a)
            mpfr_mul(&da, &alpha, &xa, MPFR_RNDD);
            mpfr_ui_div(&temp, 1, &xa, MPFR_RNDU);
            mpfr_sub(&da, &temp, &da, MPFR_RNDU);

            // compute difference (1/b - alpha b)
            mpfr_mul(&db, &alpha, &xb, MPFR_RNDD);
            mpfr_ui_div(&temp, 1, &xb, MPFR_RNDU);
            mpfr_sub(&db, &temp, &db, MPFR_RNDU);

            mpfr_max(&db, &da, &db, MPFR_RNDN);

            // compute difference (1/u - alpha u)
            mpfr_neg(&du, &alpha, MPFR_RNDN);
            mpfr_sqrt(&du, &du, MPFR_RNDD);
            mpfr_mul_si(&du, &du, 2, MPFR_RNDD);

            // compute gamma
            mpfr_add(&gamma, &db, &du, MPFR_RNDN);
            mpfr_div_si(&gamma, &gamma, 2, MPFR_RNDN);

            // compute delta
            mpfr_sub(&delta, &gamma, &du, MPFR_RNDU);
            mpfr_sub(&temp, &db, &gamma, MPFR_RNDU);
            mpfr_max(&delta, &delta, &temp, MPFR_RNDU);

            // Handle negative x.
            if (sign < 0) {
                mpfr_neg(&gamma, &gamma, MPFR_RNDN);
            }

            // compute affine approximation
            arpra_affine_1(z, x, &alpha, &gamma, &delta);
        }
    }

    // Clear vars.
    mpfr_clear(&temp);
    mpfr_clear(&xa);
    mpfr_clear(&xb);
    mpfr_clear(&da);
    mpfr_clear(&db);
    mpfr_clear(&du);
    mpfr_clear(&alpha);
    mpfr_clear(&gamma);
    mpfr_clear(&delta);
}
