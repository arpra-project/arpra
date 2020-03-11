/*
 * set.c -- Set one Arpra range with the value of another.
 *
 * Copyright 2016-2020 James Paul Turner.
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

void arpra_set (arpra_range *y, const arpra_range *x1)
{
    mpfi_t ia_range;
    mpfr_t temp1, temp2, error;
    arpra_prec prec_internal;
    arpra_uint iy, ix1;

    // Handle the y = x1 case.
    if (y == x1) return;

    // Domain violations:
    // (NaN) = (NaN)
    // (Inf) = (Inf)

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
    mpfi_init2(ia_range, y->precision);
    mpfr_init2(temp1, prec_internal + 8);
    mpfr_init2(temp2, prec_internal + 8);
    mpfr_init2(error, prec_internal);
    mpfr_set_prec(&(y->centre), prec_internal);
    mpfr_set_prec(&(y->radius), prec_internal);
    mpfr_set_zero(error, 1);
    mpfr_set_zero(&(y->radius), 1);
    arpra_helper_clear_terms(y);

    // MPFI set
    mpfi_set(ia_range, &(x1->true_range));

    // y[0] = x1[0]
    ARPRA_MPFR_RNDERR_SET(error, MPFR_RNDN, &(y->centre), &(x1->centre));

    // Allocate memory for deviation terms.
    y->symbols = malloc((x1->nTerms + 1) * sizeof(arpra_uint));
    y->deviations = malloc((x1->nTerms + 1) * sizeof(mpfr_t));

    for (iy = 0, ix1 = 0; ix1 < x1->nTerms; ix1++) {
        mpfr_init2(&(y->deviations[iy]), prec_internal);

        // y[i] = x1[i]
        y->symbols[iy] = x1->symbols[iy];
        ARPRA_MPFR_RNDERR_SET(error, MPFR_RNDN, &(y->deviations[iy]), &(x1->deviations[ix1]));

        mpfr_abs(temp1, &(y->deviations[iy]), MPFR_RNDU);
        mpfr_add(&(y->radius), &(y->radius), temp1, MPFR_RNDU);
        iy++;
    }

    // Store new deviation term.
    y->symbols[iy] = arpra_helper_next_symbol();
    y->deviations[iy] = *error;
    mpfr_add(&(y->radius), &(y->radius), &(y->deviations[iy]), MPFR_RNDU);
    y->nTerms = iy + 1;

    // Compute true_range.
    arpra_helper_compute_range(y);

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
    mpfr_clear(temp1);
    mpfr_clear(temp2);
}
