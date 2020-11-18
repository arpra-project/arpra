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
    mpfi_t ia_range_working_prec, ia_range_internal_prec, x1_range;
    mpfi_t alpha, gamma;
    mpfr_t delta;
    mpfi_t diff1, diff2, diff3;
    mpfi_srcptr diff_boundary;
    mpfi_t temp1, temp2;
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
    mpfi_init2(ia_range_working_prec, y->precision);
    mpfi_init2(ia_range_internal_prec, prec_internal);
    mpfi_init2(x1_range, x1->precision);
    mpfi_init2(alpha, prec_internal);
    mpfi_init2(gamma, prec_internal);
    mpfr_init2(delta, prec_internal);
    mpfi_init2(diff1, prec_internal);
    mpfi_init2(diff2, prec_internal);
    mpfi_init2(diff3, prec_internal);
    mpfi_init2(temp1, prec_internal);
    mpfi_init2(temp2, prec_internal);

    sign = mpfr_sgn(&(x1->true_range.left));
    if (sign < 0) {
        mpfi_neg(x1_range, &(x1->true_range));
    }
    else {
        mpfi_set(x1_range, &(x1->true_range));
    }

    mpfi_inv(ia_range_internal_prec, x1_range);

#if ARPRA_MIN_RANGE

    // compute alpha
    mpfi_set_fr(temp1, &(x1_range->right));
    mpfi_si_div(temp1, -1, temp1);
    mpfi_div_fr(alpha, temp1, &(x1_range->right));

    // compute difference (1/a - alpha a)
    mpfi_set_fr(temp1, &(ia_range_internal_prec->right));
    mpfi_mul_fr(temp2, alpha, &(x1_range->left));
    mpfi_sub(diff1, temp1, temp2);

    // compute difference (1/b - alpha b)
    mpfi_set_fr(temp1, &(ia_range_internal_prec->left));
    mpfi_mul_fr(temp2, alpha, &(x1_range->right));
    mpfi_sub(diff3, temp1, temp2);

    // compute gamma
    mpfi_add(gamma, diff1, diff3);
    mpfi_div_si(gamma, gamma, 2);

    // compute delta
    mpfi_sub(temp1, gamma, diff3);
    mpfi_sub(temp2, diff1, gamma);
    mpfr_max(delta, &(temp1->right), &(temp2->right), MPFR_RNDU);

#else

    // compute alpha
    mpfi_set_fr(temp1, &(x1_range->right));
    mpfi_si_div(temp1, -1, temp1);
    mpfi_div_fr(alpha, temp1, &(x1_range->left));

    // compute difference (1/a - alpha a)
    mpfi_set_fr(temp1, &(ia_range_internal_prec->right));
    mpfi_mul_fr(temp2, alpha, &(x1_range->left));
    mpfi_sub(diff1, temp1, temp2);

    // compute difference (1/b - alpha b)
    mpfi_set_fr(temp1, &(ia_range_internal_prec->left));
    mpfi_mul_fr(temp2, alpha, &(x1_range->right));
    mpfi_sub(diff3, temp1, temp2);

    // max of boundary diffs
    if (mpfr_greater_p(&(diff1->right), &(diff3->right))) {
        diff_boundary = diff1;
    }
    else {
        diff_boundary = diff3;
    }

    // compute difference (1/u - alpha u)
    mpfi_neg(temp1, alpha);
    mpfi_sqrt(temp1, temp1);
    mpfi_mul_si(diff2, temp1, 2);

    // compute gamma
    mpfi_add(gamma, diff_boundary, diff2);
    mpfi_div_si(gamma, gamma, 2);

    // compute delta
    mpfi_sub(temp1, gamma, diff2);
    mpfi_sub(temp2, diff_boundary, gamma);
    mpfr_max(delta, &(temp1->right), &(temp2->right), MPFR_RNDU);

#endif // ARPRA_MIN_RANGE

    if (sign < 0) {
        mpfi_neg(gamma, gamma);
    }

    // MPFI inverse
    mpfi_inv(ia_range_working_prec, &(x1->true_range));

    // compute affine approximation
    arpra_helper_affine_1(y, x1, alpha, gamma, delta);

    // Compute true_range.
    arpra_helper_compute_range(y);

    // Mix with IA range, and trim error term.
    arpra_helper_mix_trim(y, ia_range_working_prec);

    // Check for NaN and Inf.
    arpra_helper_check_result(y);

    // Clear vars.
    mpfi_clear(ia_range_working_prec);
    mpfi_clear(ia_range_internal_prec);
    mpfi_clear(x1_range);
    mpfi_clear(alpha);
    mpfi_clear(gamma);
    mpfr_clear(delta);
    mpfi_clear(diff1);
    mpfi_clear(diff2);
    mpfi_clear(diff3);
    mpfi_clear(temp1);
    mpfi_clear(temp2);
}
