/*
 * inv.c -- Compute the inverse of an Arpra range.
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
 * This affine inverse function uses a Chebyshev linear approximation.
 */

void arpra_inv (arpra_range *y, const arpra_range *x1)
{
    mpfi_t ia_range, x1_range;
    mpfr_t alpha, gamma, delta;
    mpfr_t diff1, diff2, diff3;
    mpfr_t temp1, temp2;
    arpra_prec prec_internal;
    int sign;

    // Domain violations:
    // inv(NaN) = (NaN)
    // inv(Inf) = (Inf)
    // inv(0)   = (Inf)

    // Handle domain violations.
    if (arpra_nan_p(x1)) {
        arpra_set_nan(y);
        return;
    }
    if (arpra_has_zero_p(x1)) {
        arpra_set_inf(y);
        return;
    }

    // Handle zero-width x1.
    if (mpfr_equal_p(&(x1->true_range.left), &(x1->true_range.right))) {
        arpra_mpfr_ui_fn2(mpfr_ui_div, y, 1, &(x1->true_range.left));
        return;
    }

    // Initialise vars.
    prec_internal = arpra_get_internal_precision();
    mpfi_init2(ia_range, y->precision);
    mpfi_init2(x1_range, x1->precision);
    mpfr_init2(alpha, prec_internal);
    mpfr_init2(gamma, prec_internal);
    mpfr_init2(delta, prec_internal);
    mpfr_init2(diff1, prec_internal);
    mpfr_init2(diff2, prec_internal);
    mpfr_init2(diff3, prec_internal);
    mpfr_init2(temp1, prec_internal);
    mpfr_init2(temp2, prec_internal);

    sign = mpfr_sgn(&(x1->true_range.left));
    if (sign < 0) {
        mpfr_neg(&(x1_range->left), &(x1->true_range.right), MPFR_RNDD);
        mpfr_neg(&(x1_range->right), &(x1->true_range.left), MPFR_RNDU);
    }
    else {
        mpfi_set(x1_range, &(x1->true_range));
    }

#if ARPRA_MIN_RANGE
    // compute alpha
    mpfr_si_div(alpha, -1, &(x1_range->right), MPFR_RNDN);
    mpfr_div(alpha, alpha, &(x1_range->right), MPFR_RNDN);

    // compute difference (1/a - alpha a)
    mpfr_mul(diff1, alpha, &(x1_range->left), MPFR_RNDD);
    mpfr_ui_div(temp1, 1, &(x1_range->left), MPFR_RNDU);
    mpfr_sub(diff1, temp1, diff1, MPFR_RNDU);

    // compute difference (1/b - alpha b)
    mpfr_mul(diff3, alpha, &(x1_range->right), MPFR_RNDU);
    mpfr_ui_div(temp1, 1, &(x1_range->right), MPFR_RNDD);
    mpfr_sub(diff3, temp1, diff3, MPFR_RNDD);

    // compute gamma
    mpfr_add(gamma, diff1, diff3, MPFR_RNDN);
    mpfr_div_si(gamma, gamma, 2, MPFR_RNDN);

    // compute delta
    mpfr_sub(delta, gamma, diff3, MPFR_RNDU);
    mpfr_sub(temp1, diff1, gamma, MPFR_RNDU);
    mpfr_max(delta, delta, temp1, MPFR_RNDU);

#else
    // compute alpha
    mpfr_si_div(alpha, -1, &(x1_range->right), MPFR_RNDN);
    mpfr_div(alpha, alpha, &(x1_range->left), MPFR_RNDN);

    // compute difference (1/a - alpha a)
    mpfr_mul(diff1, alpha, &(x1_range->left), MPFR_RNDD);
    mpfr_ui_div(temp1, 1, &(x1_range->left), MPFR_RNDU);
    mpfr_sub(diff1, temp1, diff1, MPFR_RNDU);

    // compute difference (1/b - alpha b)
    mpfr_mul(diff3, alpha, &(x1_range->right), MPFR_RNDD);
    mpfr_ui_div(temp1, 1, &(x1_range->right), MPFR_RNDU);
    mpfr_sub(diff3, temp1, diff3, MPFR_RNDU);

    mpfr_max(diff3, diff1, diff3, MPFR_RNDU);

    // compute difference (1/u - alpha u)
    mpfr_neg(diff2, alpha, MPFR_RNDN);
    mpfr_sqrt(diff2, diff2, MPFR_RNDD);
    mpfr_mul_si(diff2, diff2, 2, MPFR_RNDD);

    // compute gamma
    mpfr_add(gamma, diff3, diff2, MPFR_RNDN);
    mpfr_div_si(gamma, gamma, 2, MPFR_RNDN);

    // compute delta
    mpfr_sub(delta, gamma, diff2, MPFR_RNDU);
    mpfr_sub(temp1, diff3, gamma, MPFR_RNDU);
    mpfr_max(delta, delta, temp1, MPFR_RNDU);

#endif // ARPRA_MIN_RANGE

    if (sign < 0) {
        mpfr_neg(gamma, gamma, MPFR_RNDN);
    }

    // MPFI inverse
    mpfi_inv(ia_range, &(x1->true_range));

    // compute affine approximation
    arpra_affine_1(y, x1, alpha, gamma, delta);

    // Compute true_range.
    arpra_helper_compute_range(y);

    // Mix with IA range, and trim error term.
    arpra_helper_mix_trim(y, ia_range);

    // Check for NaN and Inf.
    arpra_helper_check_result(y);

    // Clear vars.
    mpfi_clear(ia_range);
    mpfi_clear(x1_range);
    mpfr_clear(alpha);
    mpfr_clear(gamma);
    mpfr_clear(delta);
    mpfr_clear(diff1);
    mpfr_clear(diff2);
    mpfr_clear(diff3);
    mpfr_clear(temp1);
    mpfr_clear(temp2);
}
