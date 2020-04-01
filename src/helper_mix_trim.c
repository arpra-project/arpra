/*
 * helper_mix_trim.c -- Mixed IA/AA and Mixed Trimmed IA/AA functions.
 *
 * Copyright 2020 James Paul Turner.
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

void arpra_helper_mix_trim (arpra_range *y, mpfi_srcptr ia_range)
{
    mpfr_t temp1, temp2;
    arpra_uint prec_internal;

    // Mixed IA/AA method.
    if (arpra_get_range_method() == ARPRA_MIXED_IAAA) {
        // Intersect AA and IA ranges.
        mpfi_intersect(&(y->true_range), &(y->true_range), ia_range);
        assert(!mpfi_is_empty(&(y->true_range)));
    }

    // Mixed trimmed IA/AA method.
    else if (arpra_get_range_method() == ARPRA_MIXED_TRIMMED_IAAA) {
        // Intersect AA and IA ranges.
        mpfi_intersect(&(y->true_range), &(y->true_range), ia_range);
        assert(!mpfi_is_empty(&(y->true_range)));

        // Initialise vars.
        prec_internal = arpra_get_internal_precision();
        mpfr_init2(temp1, prec_internal * 2);
        mpfr_init2(temp2, prec_internal * 2);

        // Trim error term if AA range fully encloses mixed IA/AA range.
        mpfr_sub(temp1, &(y->centre), &(y->radius), MPFR_RNDD);
        mpfr_add(temp2, &(y->centre), &(y->radius), MPFR_RNDU);
        if (mpfr_less_p(temp1, &(y->true_range.left))
            && mpfr_greater_p(temp2, &(y->true_range.right))) {
            mpfr_sub(temp1, &(y->true_range.left), temp1, MPFR_RNDD);
            mpfr_sub(temp2, temp2, &(y->true_range.right), MPFR_RNDD);
            mpfr_min(temp1, temp1, temp2, MPFR_RNDD);

            mpfr_sub(&(y->deviations[y->nTerms - 1]), &(y->deviations[y->nTerms - 1]), temp1, MPFR_RNDU);
            if (mpfr_sgn(&(y->deviations[y->nTerms - 1])) < 0) {
                mpfr_set_zero(&(y->deviations[y->nTerms - 1]), 1);
            }

            mpfr_sub(&(y->radius), &(y->radius), temp1, MPFR_RNDU);
            if (mpfr_sgn(&(y->radius)) < 0) {
                mpfr_set_zero(&(y->radius), 1);
            }
        }

        // Clear vars.
        mpfr_clear(temp1);
        mpfr_clear(temp2);
    }
}
