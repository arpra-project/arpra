/*
 * log.c -- Compute the natural logarithm of an Arpra range.
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
 * This affine natural logarithm function uses a Chebyshev linear approximation.
 */

void arpra_log (arpra_range *y, const arpra_range *x1)
{
    mpfi_t ia_range_working_prec, ia_range_internal_prec;
    mpfi_t alpha, gamma;
    mpfr_t delta;
    mpfi_t diff1, diff2, diff3;
    mpfi_srcptr diff_lo, diff_hi;
    mpfi_t temp1, temp2;
    arpra_prec prec_internal;

    // Domain violations:
    // log(NaN)   = (NaN)
    // log(Inf)   = (NaN)
    // log(R < 0) = (NaN)

    // Handle domain violations.
    if (arpra_nan_p(x1) || arpra_has_neg_p(x1)) {
        arpra_set_nan(y);
        return;
    }

    // Handle zero-width x1.
    if (mpfr_equal_p(&(x1->true_range.left), &(x1->true_range.right))) {
        arpra_mpfr_fn1(mpfr_log, y, &(x1->true_range.left));
        return;
    }

    // Initialise vars.
    prec_internal = arpra_get_internal_precision();
    mpfi_init2(ia_range_working_prec, y->precision);
    mpfi_init2(ia_range_internal_prec, prec_internal);
    mpfi_init2(alpha, prec_internal);
    mpfi_init2(gamma, prec_internal);
    mpfr_init2(delta, prec_internal);
    mpfi_init2(diff1, prec_internal);
    mpfi_init2(diff2, prec_internal);
    mpfi_init2(diff3, prec_internal);
    mpfi_init2(temp1, prec_internal);
    mpfi_init2(temp2, prec_internal);

    mpfi_log(ia_range_internal_prec, &(x1->true_range));

    // compute alpha
    mpfi_set_fr(temp1, &(ia_range_internal_prec->left));
    mpfi_set_fr(temp2, &(ia_range_internal_prec->right));
    mpfi_sub(alpha, temp2, temp1);
    mpfi_set_fr(temp1, &(x1->true_range.left));
    mpfi_set_fr(temp2, &(x1->true_range.right));
    mpfi_sub(temp1, temp2, temp1);
    mpfi_div(alpha, alpha, temp1);

    // compute difference (log(a) - alpha a)
    mpfi_set_fr(temp1, &(ia_range_internal_prec->left));
    mpfi_mul_fr(temp2, alpha, &(x1->true_range.left));
    mpfi_sub(diff1, temp1, temp2);

    // compute difference (log(b) - alpha b)
    mpfi_set_fr(temp1, &(ia_range_internal_prec->right));
    mpfi_mul_fr(temp2, alpha, &(x1->true_range.right));
    mpfi_sub(diff3, temp1, temp2);

    // compute difference (log(u) - alpha u)
    mpfi_si_div(diff2, 1, alpha);
    mpfi_log(diff2, diff2);
    mpfi_sub_si(diff2, diff2, 1);

    // min and max difference
    diff_lo = mpfr_less_p(&(diff1->left), &(diff3->left)) ? diff1 : diff3;
    diff_hi = diff2;

    // compute gamma
    mpfi_add(gamma, diff_lo, diff_hi);
    mpfi_div_si(gamma, gamma, 2);

    // compute delta
    mpfi_sub(temp1, gamma, diff_lo);
    mpfi_sub(temp2, diff_hi, gamma);
    mpfr_max(delta, &(temp1->right), &(temp2->right), MPFR_RNDU);

    // MPFI natural logarithm
    mpfi_log(ia_range_working_prec, &(x1->true_range));

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
    mpfi_clear(alpha);
    mpfi_clear(gamma);
    mpfr_clear(delta);
    mpfi_clear(diff1);
    mpfi_clear(diff2);
    mpfi_clear(diff3);
    mpfi_clear(temp1);
    mpfi_clear(temp2);
}
