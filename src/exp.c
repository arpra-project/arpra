/*
 * exp.c -- Compute the exponent of an Arpra range.
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
 * This affine exponential function uses a Chebyshev linear approximation.
 */

void arpra_exp (arpra_range *y, const arpra_range *x1)
{
    mpfi_t ia_range;
    mpfr_t alpha, gamma, delta;
    mpfr_t diff1, diff2, diff3;
    mpfr_t temp1, temp2;
    arpra_prec prec_internal;

    // Domain violations:
    // exp(NaN) = (NaN)
    // exp(Inf) = (Inf)

    // Handle domain violations.
    if (arpra_nan_p(x1)) {
        arpra_set_nan(y);
        return;
    }
    if (arpra_inf_p(x1)) {
        arpra_set_inf(y);
        return;
    }

    // Handle zero-width x1.
    if (mpfr_equal_p(&(x1->true_range.left), &(x1->true_range.right))) {
        arpra_mpfr_f1(mpfr_exp, y, &(x1->true_range.left));
        return;
    }

    // Initialise vars.
    prec_internal = arpra_get_internal_precision();
    mpfi_init2(ia_range, y->precision);
    mpfr_init2(alpha, prec_internal);
    mpfr_init2(gamma, prec_internal);
    mpfr_init2(delta, prec_internal);
    mpfr_init2(diff1, prec_internal);
    mpfr_init2(diff2, prec_internal);
    mpfr_init2(diff3, prec_internal);
    mpfr_init2(temp1, prec_internal);
    mpfr_init2(temp2, prec_internal);

#if ARPRA_MIN_RANGE
    // compute alpha
    mpfr_exp(alpha, &(x1->true_range.left), MPFR_RNDN);

    // compute difference (exp(a) - alpha a)
    mpfr_mul(diff1, alpha, &(x1->true_range.left), MPFR_RNDU);
    mpfr_exp(temp1, &(x1->true_range.left), MPFR_RNDD);
    mpfr_sub(diff1, temp1, diff1, MPFR_RNDD);

    // compute difference (exp(b) - alpha b)
    mpfr_mul(diff3, alpha, &(x1->true_range.right), MPFR_RNDD);
    mpfr_exp(temp1, &(x1->true_range.right), MPFR_RNDU);
    mpfr_sub(diff3, temp1, diff3, MPFR_RNDU);

    // compute gamma
    mpfr_add(gamma, diff1, diff3, MPFR_RNDN);
    mpfr_div_si(gamma, gamma, 2, MPFR_RNDN);

    // compute delta
    mpfr_sub(delta, gamma, diff1, MPFR_RNDU);
    mpfr_sub(temp1, diff3, gamma, MPFR_RNDU);
    mpfr_max(delta, delta, temp1, MPFR_RNDU);

#else
    // compute alpha
    mpfr_exp(alpha, &(x1->true_range.right), MPFR_RNDN);
    mpfr_exp(temp1, &(x1->true_range.left), MPFR_RNDN);
    mpfr_sub(alpha, alpha, temp1, MPFR_RNDN);
    mpfr_sub(temp1, &(x1->true_range.right), &(x1->true_range.left), MPFR_RNDN);
    mpfr_div(alpha, alpha, temp1, MPFR_RNDN);

    // compute difference (exp(a) - alpha a)
    mpfr_mul(diff1, alpha, &(x1->true_range.left), MPFR_RNDD);
    mpfr_exp(temp1, &(x1->true_range.left), MPFR_RNDU);
    mpfr_sub(diff1, temp1, diff1, MPFR_RNDU);

    // compute difference (exp(b) - alpha b)
    mpfr_mul(diff3, alpha, &(x1->true_range.right), MPFR_RNDD);
    mpfr_exp(temp1, &(x1->true_range.right), MPFR_RNDU);
    mpfr_sub(diff3, temp1, diff3, MPFR_RNDU);

    mpfr_max(diff3, diff1, diff3, MPFR_RNDU);

    // compute difference (exp(u) - alpha u)
    mpfr_log(diff2, alpha, MPFR_RNDU);
    mpfr_sub_si(diff2, diff2, 1, MPFR_RNDU);
    mpfr_mul(diff2, alpha, diff2, MPFR_RNDU);
    mpfr_neg(diff2, diff2, MPFR_RNDD);

    // compute gamma
    mpfr_add(gamma, diff3, diff2, MPFR_RNDN);
    mpfr_div_si(gamma, gamma, 2, MPFR_RNDN);

    // compute delta
    mpfr_sub(delta, gamma, diff2, MPFR_RNDU);
    mpfr_sub(temp1, diff3, gamma, MPFR_RNDU);
    mpfr_max(delta, delta, temp1, MPFR_RNDU);

#endif // ARPRA_MIN_RANGE

    // MPFI exponential
    mpfi_exp(ia_range, &(x1->true_range));

    // compute affine approximation
    arpra_affine_1(y, x1, alpha, gamma, delta);

    // Compute true_range, and add rounding error.
    arpra_helper_range_rounded(y);

#ifdef ARPRA_MIXED_IAAA
    // Intersect AA and IA ranges.
    mpfi_intersect(&(y->true_range), &(y->true_range), ia_range);

#ifdef ARPRA_MIXED_TRIMMED_IAAA
    // Trim error term if AA range fully encloses mixed IA/AA range.
    mpfr_sub(temp1, &(y->centre), &(y->radius), MPFR_RNDD);
    mpfr_add(temp2, &(y->centre), &(y->radius), MPFR_RNDU);
    if (mpfr_less_p(temp1, &(y->true_range.left))
        && mpfr_greater_p(temp2, &(y->true_range.right))) {
        mpfr_sub(temp1, &(y->true_range.left), temp1, MPFR_RNDD);
        mpfr_sub(temp2, temp2, &(y->true_range.right), MPFR_RNDD);
        mpfr_min(temp1, temp1, temp2, MPFR_RNDD);
        if (mpfr_greater_p(temp1, &(y->deviations[y->nTerms - 1]))) {
            mpfr_sub(&(y->radius), &(y->radius), &(y->deviations[y->nTerms - 1]), MPFR_RNDU);
            mpfr_set_zero(&(y->deviations[y->nTerms - 1]), 1);
        }
        else {
            mpfr_sub(&(y->radius), &(y->radius), temp1, MPFR_RNDU);
            mpfr_sub(&(y->deviations[y->nTerms - 1]), &(y->deviations[y->nTerms - 1]), temp1, MPFR_RNDU);
        }
    }
#endif // ARPRA_MIXED_TRIMMED_IAAA
#endif // ARPRA_MIXED_IAAA

    // Check for NaN and Inf.
    arpra_helper_check_result(y);

    // Clear vars.
    mpfi_clear(ia_range);
    mpfr_clear(alpha);
    mpfr_clear(gamma);
    mpfr_clear(delta);
    mpfr_clear(diff1);
    mpfr_clear(diff2);
    mpfr_clear(diff3);
    mpfr_clear(temp1);
    mpfr_clear(temp2);
}
