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
    arpra_mpfr alpha, gamma, delta;
    arpra_mpfr temp1, temp2, da, db, du;
    arpra_mpfi ia_range, x_range;
    arpra_prec prec_internal;
    arpra_int sign;

    // Initialise vars.
    prec_internal = arpra_get_internal_precision();
    mpfr_init2(&alpha, prec_internal);
    mpfr_init2(&gamma, prec_internal);
    mpfr_init2(&delta, prec_internal);
    mpfr_init2(&temp1, prec_internal);
    mpfr_init2(&temp2, prec_internal);
    mpfr_init2(&da, prec_internal);
    mpfr_init2(&db, prec_internal);
    mpfr_init2(&du, prec_internal);
    mpfi_init2(&ia_range, z->precision);
    mpfi_init2(&x_range, x->precision);

    // Handle x with zero radius.
    if (mpfr_zero_p(&(x->radius))) {
        if (mpfr_ui_div(&temp1, 1, &(x->centre), MPFR_RNDN)) {
            arpra_helper_error_half_ulp(&delta, &temp1);
            arpra_set_mpfr_rad(z, &temp1, &delta);
        }
        else {
            arpra_set_mpfr(z, &temp1);
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
            // MPFI inverse
            mpfi_inv(&ia_range, &(x->true_range));

            sign = mpfr_sgn(&(x->true_range.left));
            if (sign < 0) {
                mpfr_neg(&(x_range.left), &(x->true_range.right), MPFR_RNDD);
                mpfr_neg(&(x_range.right), &(x->true_range.left), MPFR_RNDU);
            }
            else {
                mpfi_set(&x_range, &(x->true_range));
            }

#if ARPRA_MIN_RANGE
            // compute alpha
            mpfr_si_div(&alpha, -1, &(x_range.right), MPFR_RNDN);
            mpfr_div(&alpha, &alpha, &(x_range.right), MPFR_RNDN);

            // compute difference (1/a - alpha a)
            mpfr_mul(&da, &alpha, &(x_range.left), MPFR_RNDD);
            mpfr_ui_div(&temp1, 1, &(x_range.left), MPFR_RNDU);
            mpfr_sub(&da, &temp1, &da, MPFR_RNDU);

            // compute difference (1/b - alpha b)
            mpfr_mul(&db, &alpha, &(x_range.right), MPFR_RNDU);
            mpfr_ui_div(&temp1, 1, &(x_range.right), MPFR_RNDD);
            mpfr_sub(&db, &temp1, &db, MPFR_RNDD);

            // compute gamma
            mpfr_add(&gamma, &da, &db, MPFR_RNDN);
            mpfr_div_si(&gamma, &gamma, 2, MPFR_RNDN);

            // compute delta
            mpfr_sub(&delta, &gamma, &db, MPFR_RNDU);
            mpfr_sub(&temp1, &da, &gamma, MPFR_RNDU);
            mpfr_max(&delta, &delta, &temp1, MPFR_RNDU);

#else
            // compute alpha
            mpfr_si_div(&alpha, -1, &(x_range.right), MPFR_RNDN);
            mpfr_div(&alpha, &alpha, &(x_range.left), MPFR_RNDN);

            // compute difference (1/a - alpha a)
            mpfr_mul(&da, &alpha, &(x_range.left), MPFR_RNDD);
            mpfr_ui_div(&temp1, 1, &(x_range.left), MPFR_RNDU);
            mpfr_sub(&da, &temp1, &da, MPFR_RNDU);

            // compute difference (1/b - alpha b)
            mpfr_mul(&db, &alpha, &(x_range.right), MPFR_RNDD);
            mpfr_ui_div(&temp1, 1, &(x_range.right), MPFR_RNDU);
            mpfr_sub(&db, &temp1, &db, MPFR_RNDU);

            mpfr_max(&db, &da, &db, MPFR_RNDU);

            // compute difference (1/u - alpha u)
            mpfr_neg(&du, &alpha, MPFR_RNDN);
            mpfr_sqrt(&du, &du, MPFR_RNDD);
            mpfr_mul_si(&du, &du, 2, MPFR_RNDD);

            // compute gamma
            mpfr_add(&gamma, &db, &du, MPFR_RNDN);
            mpfr_div_si(&gamma, &gamma, 2, MPFR_RNDN);

            // compute delta
            mpfr_sub(&delta, &gamma, &du, MPFR_RNDU);
            mpfr_sub(&temp1, &db, &gamma, MPFR_RNDU);
            mpfr_max(&delta, &delta, &temp1, MPFR_RNDU);

#endif // ARPRA_MIN_RANGE

            if (sign < 0) {
                mpfr_neg(&gamma, &gamma, MPFR_RNDN);
            }

            // compute affine approximation
            arpra_affine_1(z, x, &alpha, &gamma, &delta);

#ifdef ARPRA_MIXED_IAAA
            // Intersect AA and IA ranges.
            mpfi_intersect(&(z->true_range), &(z->true_range), &ia_range);

#ifdef ARPRA_MIXED_TRIMMED_IAAA
            // Trim error term if AA range fully encloses mixed IA/AA range.
            mpfr_sub(&temp1, &(z->centre), &(z->radius), MPFR_RNDD);
            mpfr_add(&temp2, &(z->centre), &(z->radius), MPFR_RNDU);
            if (mpfr_less_p(&temp1, &(z->true_range.left))
                && mpfr_greater_p(&temp2, &(z->true_range.right))) {
                mpfr_sub(&temp1, &(z->true_range.left), &temp1, MPFR_RNDD);
                mpfr_sub(&temp2, &temp2, &(z->true_range.right), MPFR_RNDD);
                mpfr_min(&temp1, &temp1, &temp2, MPFR_RNDD);
                if (mpfr_greater_p(&temp1, &(z->deviations[z->nTerms - 1]))) {
                    mpfr_sub(&(z->radius), &(z->radius), &(z->deviations[z->nTerms - 1]), MPFR_RNDU);
                    mpfr_set_zero(&(z->deviations[z->nTerms - 1]), 1);
                }
                else {
                    mpfr_sub(&(z->radius), &(z->radius), &temp1, MPFR_RNDU);
                    mpfr_sub(&(z->deviations[z->nTerms - 1]), &(z->deviations[z->nTerms - 1]), &temp1, MPFR_RNDU);
                }
            }
#endif // ARPRA_MIXED_TRIMMED_IAAA
#endif // ARPRA_MIXED_IAAA
        }
    }

    // Clear vars.
    mpfr_clear(&alpha);
    mpfr_clear(&gamma);
    mpfr_clear(&delta);
    mpfr_clear(&temp1);
    mpfr_clear(&temp2);
    mpfr_clear(&da);
    mpfr_clear(&db);
    mpfr_clear(&du);
    mpfi_clear(&ia_range);
    mpfi_clear(&x_range);
}
