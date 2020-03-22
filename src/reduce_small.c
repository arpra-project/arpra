/*
 * reduce_small.c -- Reduce deviation terms less than a given threshold.
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

void arpra_reduce_small_abs (arpra_range *y, const arpra_range *x1, mpfr_srcptr abs_threshold)
{
    mpfr_t error;
    mpfr_ptr sum_x, *sum_x_ptr;
    arpra_range yy;
    arpra_prec prec_internal;
    arpra_uint i_y, i_x1;

    // Handle trivial cases.
    if (mpfr_sgn(abs_threshold) < 0) {
        arpra_set(y, x1);
        return;
    }

    // Domain violations:
    // reduce(NaN) = (NaN)
    // reduce(Inf) = (Inf)

    // Handle domain violations.
    if (arpra_nan_p(x1)) {
        arpra_set_nan(y);
        return;
    }
    if (arpra_inf_p(x1)) {
        arpra_set_inf(y);
        return;
    }

    // Initialise vars.
    prec_internal = arpra_get_internal_precision();
    mpfr_init2(error, prec_internal);
    arpra_init2(&yy, y->precision);
    sum_x = malloc((x1->nTerms + 1) * sizeof(mpfr_t));
    sum_x_ptr = malloc((x1->nTerms + 1) * sizeof(mpfr_ptr));
    mpfr_set_zero(error, 1);

    // y[0] = x1[0]
    ARPRA_MPFR_RNDERR_SET(error, MPFR_RNDN, &(yy.centre), &(x1->centre));

    // Allocate memory for deviation terms.
    yy.symbols = malloc((x1->nTerms + 1) * sizeof(arpra_uint));
    yy.deviations = malloc((x1->nTerms + 1) * sizeof(mpfr_t));

    for (i_y = 0, i_x1 = 0; i_x1 < x1->nTerms; i_x1++) {
        if (mpfr_cmpabs(&(x1->deviations[i_x1]), abs_threshold) > 0) {
            mpfr_init2(&(yy.deviations[i_y]), prec_internal);

            // y[i] = x1[i]
            yy.symbols[i_y] = x1->symbols[i_x1];
            ARPRA_MPFR_RNDERR_SET(error, MPFR_RNDN, &(yy.deviations[i_y]), &(x1->deviations[i_x1]));

            i_y++;
        }
        else {
            // Include small deviation term in abs sum.
            sum_x[i_x1 - i_y] = x1->deviations[i_x1];
            sum_x[i_x1 - i_y]._mpfr_sign = 1;
            sum_x_ptr[i_x1 - i_y] = &(sum_x[i_x1 - i_y]);
        }
    }

    // Abs sum the small deviation terms.
    sum_x_ptr[i_x1 - i_y] = error;
    mpfr_sum(error, sum_x_ptr, (i_x1 - i_y + 1), MPFR_RNDU);

    // Store new deviation term.
    yy.symbols[i_y] = arpra_helper_next_symbol();
    yy.deviations[i_y] = *error;
    yy.nTerms = i_y + 1;

    // Compute true_range.
    arpra_helper_compute_range(&yy);

    // Mix with IA range, and trim error term.
    arpra_helper_mix_trim(&yy, &(x1->true_range));

    // Check for NaN and Inf.
    arpra_helper_check_result(&yy);

    // Clear vars.
    arpra_clear(y);
    *y = yy;
    free(sum_x);
    free(sum_x_ptr);
}

void arpra_reduce_small_rel (arpra_range *y, const arpra_range *x1, mpfr_srcptr rel_threshold)
{
    mpfr_t abs_threshold;
    arpra_prec prec_internal;

    // Initialise vars.
    prec_internal = arpra_get_internal_precision();
    mpfr_init2(abs_threshold, prec_internal);

    // Convert to absolute threshold.
    mpfr_mul(abs_threshold, &(x1->radius), rel_threshold, MPFR_RNDU);
    arpra_reduce_small_abs(y, x1, abs_threshold);

    // Clear vars.
    mpfr_clear(abs_threshold);
}
