/*
 * log.c -- Compute the natural logarithm of an Arpra range.
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
 * This affine natural logarithm function uses a Chebyshev linear approximation.
 */

void arpra_log (arpra_range *z, const arpra_range *x)
{
    arpra_mpfr alpha, gamma, delta;
    arpra_mpfr temp1, temp2, da, db, du;
    arpra_mpfi ia_range, x_range;
    arpra_prec prec_internal;

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
        if (mpfr_log(&temp1, &(x->centre), MPFR_RNDN)) {
            arpra_helper_error_half_ulp(&delta, &temp1);
            arpra_set_mpfr_rad(z, &temp1, &delta);
        }
        else {
            arpra_set_mpfr(z, &temp1);
        }
    }
    else {
        // Handle domain violations.
        if (arpra_nan_p(x) || arpra_has_neg_p(x)) {
            arpra_set_nan(z);
        }
        else if (arpra_has_zero_p(x)) {
            arpra_set_inf(z);
        }

        // Domain is OK.
        else {
            // MPFI natural logarithm
            mpfi_log(&ia_range, &(x->true_range));

            mpfi_set(&x_range, &(x->true_range));

            // compute alpha
            mpfr_log(&alpha, &(x_range.right), MPFR_RNDN);
            mpfr_log(&temp1, &(x_range.left), MPFR_RNDN);
            mpfr_sub(&alpha, &alpha, &temp1, MPFR_RNDN);
            mpfr_sub(&temp1, &(x_range.right), &(x_range.left), MPFR_RNDN);
            mpfr_div(&alpha, &alpha, &temp1, MPFR_RNDN);

            // compute difference (log(a) - alpha a)
            mpfr_mul(&da, &alpha, &(x_range.left), MPFR_RNDU);
            mpfr_log(&temp1, &(x_range.left), MPFR_RNDD);
            mpfr_sub(&da, &temp1, &da, MPFR_RNDD);

            // compute difference (log(b) - alpha b)
            mpfr_mul(&db, &alpha, &(x_range.right), MPFR_RNDU);
            mpfr_log(&temp1, &(x_range.right), MPFR_RNDD);
            mpfr_sub(&db, &temp1, &db, MPFR_RNDD);

            mpfr_min(&da, &da, &db, MPFR_RNDD);

            // compute difference (log(u) - alpha u)
            mpfr_ui_div(&du, 1, &alpha, MPFR_RNDU);
            mpfr_log(&du, &du, MPFR_RNDU);
            mpfr_sub_si(&du, &du, 1, MPFR_RNDU);

            // compute gamma
            mpfr_add(&gamma, &da, &du, MPFR_RNDN);
            mpfr_div_si(&gamma, &gamma, 2, MPFR_RNDN);

            // compute delta
            mpfr_sub(&delta, &du, &gamma, MPFR_RNDU);
            mpfr_sub(&temp1, &gamma, &da, MPFR_RNDU);
            mpfr_max(&delta, &delta, &temp1, MPFR_RNDU);

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
                    mpfr_set_ui(&(z->deviations[z->nTerms - 1]), 0, MPFR_RNDZ);
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
