/*
 * set_mpfi.c -- Set an Arpra range with an MPFI interval.
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

void arpra_set_mpfi (arpra_range *y, mpfi_srcptr x1)
{
    mpfr_t temp1, temp2;
    arpra_prec prec_internal;

    // Domain violations:
    // (NaN) = (NaN)
    // (Inf) = (Inf)

    // Handle domain violations.
    if (mpfi_nan_p(x1)) {
        arpra_set_nan(y);
        return;
    }
    if (mpfi_inf_p(x1)) {
        arpra_set_inf(y);
        return;
    }

    // Initialise vars.
    prec_internal = arpra_get_internal_precision();
    mpfr_init2(temp1, prec_internal);
    mpfr_init2(temp2, prec_internal);
    mpfr_set_prec(&(y->centre), prec_internal);
    mpfr_set_prec(&(y->radius), prec_internal);
    arpra_helper_clear_terms(y);

    // MPFI set
    mpfi_set(&(y->true_range), x1);

    // y[0] = (x1[lo] + x1[hi]) / 2
    mpfi_mid(&(y->centre), &(y->true_range));

    // Allocate memory for deviation terms.
    y->symbols = malloc(sizeof(arpra_uint));
    y->deviations = malloc(sizeof(mpfr_t));

    // rad(y) = max{(y[0] - x1[lo]), (x1[hi] - y[0])}
    mpfr_sub(temp1, &(y->centre), &(y->true_range.left), MPFR_RNDU);
    mpfr_sub(temp2, &(y->true_range.right), &(y->centre), MPFR_RNDU);
    mpfr_max(&(y->radius), temp1, temp2, MPFR_RNDU);

    // Store new deviation term.
    y->symbols[0] = arpra_helper_next_symbol();
    mpfr_init2(&(y->deviations[0]), prec_internal);
    mpfr_set(&(y->deviations[0]), &(y->radius), MPFR_RNDU);
    y->nTerms = 1;

    // Check for NaN and Inf.
    arpra_helper_check_result(y);

    // Clear vars.
    mpfr_clear(temp1);
    mpfr_clear(temp2);
}
