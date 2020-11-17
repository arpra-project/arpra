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
    mpfi_t ia_range_working_prec, ia_range_internal_prec;
    mpfi_t alpha, gamma;
    mpfr_t delta;
    mpfi_t diff1, diff2, diff3;
    mpfi_srcptr diff_boundary;
    mpfi_t temp1, temp2;
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
        arpra_mpfr_fn1(mpfr_exp, y, &(x1->true_range.left));
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

//#if ARPRA_MIN_RANGE
    /* // compute alpha */
    /* mpfr_exp(alpha, &(x1->true_range.left), MPFR_RNDN); */

    /* // compute difference (exp(a) - alpha a) */
    /* mpfr_mul(diff1, alpha, &(x1->true_range.left), MPFR_RNDU); */
    /* mpfr_exp(temp1, &(x1->true_range.left), MPFR_RNDD); */
    /* mpfr_sub(diff1, temp1, diff1, MPFR_RNDD); */

    /* // compute difference (exp(b) - alpha b) */
    /* mpfr_mul(diff3, alpha, &(x1->true_range.right), MPFR_RNDD); */
    /* mpfr_exp(temp1, &(x1->true_range.right), MPFR_RNDU); */
    /* mpfr_sub(diff3, temp1, diff3, MPFR_RNDU); */

    /* // compute gamma */
    /* mpfr_add(gamma, diff1, diff3, MPFR_RNDN); */
    /* mpfr_div_si(gamma, gamma, 2, MPFR_RNDN); */

    /* // compute delta */
    /* mpfr_sub(delta, gamma, diff1, MPFR_RNDU); */
    /* mpfr_sub(temp1, diff3, gamma, MPFR_RNDU); */
    /* mpfr_max(delta, delta, temp1, MPFR_RNDU); */

//#else

    mpfi_exp(ia_range_internal_prec, &(x1->true_range));

    // compute alpha
    mpfi_set_fr(temp1, &(ia_range_internal_prec->left));
    mpfi_set_fr(temp2, &(ia_range_internal_prec->right));
    mpfi_sub(alpha, temp2, temp1);
    mpfi_set_fr(temp1, &(x1->true_range.left));
    mpfi_set_fr(temp2, &(x1->true_range.right));
    mpfi_sub(temp1, temp2, temp1);
    mpfi_div(alpha, alpha, temp1);

    // compute difference (exp(a) - alpha a)
    mpfi_set_fr(temp1, &(ia_range_internal_prec->left));
    mpfi_mul_fr(temp2, alpha, &(x1->true_range.left));
    mpfi_sub(diff1, temp1, temp2);

    // compute difference (exp(b) - alpha b)
    mpfi_set_fr(temp1, &(ia_range_internal_prec->right));
    mpfi_mul_fr(temp2, alpha, &(x1->true_range.right));
    mpfi_sub(diff3, temp1, temp2);

    // max of boundary diffs
    if (mpfr_greater_p(&(diff1->right), &(diff3->right))) {
        diff_boundary = diff1;
    }
    else {
        diff_boundary = diff3;
    }

    // compute difference (exp(u) - alpha u)
    mpfi_log(diff2, alpha);
    mpfi_sub_si(diff2, diff2, 1);
    mpfi_mul(diff2, alpha, diff2);
    mpfi_neg(diff2, diff2);

    // compute gamma
    mpfi_add(gamma, diff_boundary, diff2);
    mpfi_div_si(gamma, gamma, 2);

    // compute delta
    mpfi_sub(temp1, gamma, diff2);
    mpfi_sub(temp2, diff_boundary, gamma);
    mpfr_max(delta, &(temp1->right), &(temp2->right), MPFR_RNDU);

//#endif // ARPRA_MIN_RANGE

    // MPFI exponential
    mpfi_exp(ia_range_working_prec, &(x1->true_range));

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
