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
    mpfi_t ia_range;
    mpfr_t alpha, gamma, delta;
    mpfr_t diff1, diff2, diff3;
    mpfr_t temp1, temp2;
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
    mpfi_init2(ia_range, y->precision);
    mpfr_init2(alpha, prec_internal);
    mpfr_init2(gamma, prec_internal);
    mpfr_init2(delta, prec_internal);
    mpfr_init2(diff1, prec_internal);
    mpfr_init2(diff2, prec_internal);
    mpfr_init2(diff3, prec_internal);
    mpfr_init2(temp1, prec_internal);
    mpfr_init2(temp2, prec_internal);

    // compute alpha
    mpfr_log(alpha, &(x1->true_range.right), MPFR_RNDN);
    mpfr_log(temp1, &(x1->true_range.left), MPFR_RNDN);
    mpfr_sub(alpha, alpha, temp1, MPFR_RNDN);
    mpfr_sub(temp1, &(x1->true_range.right), &(x1->true_range.left), MPFR_RNDN);
    mpfr_div(alpha, alpha, temp1, MPFR_RNDN);

    // compute difference (log(a) - alpha a)
    mpfr_mul(diff1, alpha, &(x1->true_range.left), MPFR_RNDU);
    mpfr_log(temp1, &(x1->true_range.left), MPFR_RNDD);
    mpfr_sub(diff1, temp1, diff1, MPFR_RNDD);

    // compute difference (log(b) - alpha b)
    mpfr_mul(diff3, alpha, &(x1->true_range.right), MPFR_RNDU);
    mpfr_log(temp1, &(x1->true_range.right), MPFR_RNDD);
    mpfr_sub(diff3, temp1, diff3, MPFR_RNDD);

    mpfr_min(diff1, diff1, diff3, MPFR_RNDD);

    // compute difference (log(u) - alpha u)
    mpfr_ui_div(diff2, 1, alpha, MPFR_RNDU);
    mpfr_log(diff2, diff2, MPFR_RNDU);
    mpfr_sub_si(diff2, diff2, 1, MPFR_RNDU);

    // compute gamma
    mpfr_add(gamma, diff1, diff2, MPFR_RNDN);
    mpfr_div_si(gamma, gamma, 2, MPFR_RNDN);

    // compute delta
    mpfr_sub(delta, diff2, gamma, MPFR_RNDU);
    mpfr_sub(temp1, gamma, diff1, MPFR_RNDU);
    mpfr_max(delta, delta, temp1, MPFR_RNDU);

    // MPFI natural logarithm
    mpfi_log(ia_range, &(x1->true_range));

    // compute affine approximation
    arpra_helper_affine_1(y, x1, alpha, gamma, delta);

    // Compute true_range.
    arpra_helper_compute_range(y);

    // Mix with IA range, and trim error term.
    arpra_helper_mix_trim(y, ia_range);

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
