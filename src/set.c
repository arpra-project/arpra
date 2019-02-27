/*
 * set.c -- Set one Arpra range with the values of another.
 *
 * Copyright 2016-2018 James Paul Turner.
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

void arpra_set (arpra_range *z, const arpra_range *x)
{
    arpra_uint xTerm, zTerm;
    arpra_mpfr temp1, temp2, error;
    arpra_mpfi ia_range;
    arpra_prec prec_internal;

    // Handle trivial cases.
    if (z == x) return;

    // Domain violations:
    // NaN  =  NaN
    // Inf  =  Inf

    // Handle domain violations.
    if (arpra_nan_p(x)) {
        arpra_set_nan(z);
        return;
    }
    if (arpra_inf_p(x)) {
        arpra_set_inf(z);
        return;
    }

    // Initialise vars.
    prec_internal = arpra_get_internal_precision();
    mpfr_init2(&temp1, prec_internal);
    mpfr_init2(&temp2, prec_internal);
    mpfr_init2(&error, prec_internal);
    mpfi_init2(&ia_range, z->precision);
    mpfr_set_prec(&(z->centre), prec_internal);
    mpfr_set_prec(&(z->radius), prec_internal);
    mpfr_set_ui(&error, 0, MPFR_RNDU);
    mpfr_set_ui(&(z->radius), 0, MPFR_RNDU);

    // MPFI set
    mpfi_set(&ia_range, &(x->true_range));

    // z_0 = x_0
    if (mpfr_set(&(z->centre), &(x->centre), MPFR_RNDN)) {
        arpra_helper_error_half_ulp(&temp1, &(z->centre));
        mpfr_add(&error, &error, &temp1, MPFR_RNDU);
    }

    // Replace existing deviation term memory.
    arpra_clear_terms(z);
    z->nTerms = x->nTerms + 1;
    z->symbols = malloc(z->nTerms * sizeof(arpra_uint));
    z->deviations = malloc(z->nTerms * sizeof(arpra_mpfr));

    // Copy deviation terms over.
    for (xTerm = 0, zTerm = 0; xTerm < x->nTerms; xTerm++) {
        z->symbols[zTerm] = x->symbols[zTerm];
        mpfr_init2(&(z->deviations[zTerm]), prec_internal);

        // z_i = x_i
        if (mpfr_set(&(z->deviations[zTerm]), &(x->deviations[xTerm]), MPFR_RNDN)) {
            arpra_helper_error_half_ulp(&temp1, &(z->deviations[zTerm]));
            mpfr_add(&error, &error, &temp1, MPFR_RNDU);
        }

        // Store nonzero deviation terms.
        if (mpfr_zero_p(&(z->deviations[zTerm]))) {
            mpfr_clear(&(z->deviations[zTerm]));
        }
        else {
            mpfr_abs(&temp1, &(z->deviations[zTerm]), MPFR_RNDU);
            mpfr_add(&(z->radius), &(z->radius), &temp1, MPFR_RNDU);
            zTerm++;
        }
    }

    // Round range to target precision.
    mpfr_sub(&temp1, &(z->centre), &(z->radius), MPFR_RNDD);
    mpfr_sub(&temp1, &temp1, &error, MPFR_RNDD);
    mpfr_set(&(z->true_range.left), &temp1, MPFR_RNDD);
    mpfr_sub(&temp1, &temp1, &(z->true_range.left), MPFR_RNDU);
    mpfr_add(&temp2, &(z->centre), &(z->radius), MPFR_RNDU);
    mpfr_add(&temp2, &temp2, &error, MPFR_RNDU);
    mpfr_set(&(z->true_range.right), &temp2, MPFR_RNDU);
    mpfr_sub(&temp2, &(z->true_range.right), &temp2, MPFR_RNDU);
    mpfr_max(&temp1, &temp1, &temp2, MPFR_RNDU);
    mpfr_add(&error, &error, &temp1, MPFR_RNDU);

    // Store numerical error term.
    z->symbols[zTerm] = arpra_next_symbol();
    z->deviations[zTerm] = error;
    mpfr_add(&(z->radius), &(z->radius), &(z->deviations[zTerm]), MPFR_RNDU);
    z->nTerms = zTerm + 1;

#ifdef ARPRA_MIXED_IAAA
#ifdef ARPRA_MIXED_TRIMMED_IAAA
    // Trim error term if Arpra range fully contains IA range.
    if (mpfr_less_p(&(z->true_range.left), &(ia_range.left))
        && mpfr_greater_p(&(z->true_range.right), &(ia_range.right))) {
        mpfr_sub(&temp1, &(ia_range.left), &(z->true_range.left), MPFR_RNDD);
        mpfr_sub(&temp2, &(z->true_range.right), &(ia_range.right), MPFR_RNDD);
        mpfr_min(&temp1, &temp1, &temp2, MPFR_RNDD);
        mpfr_sub(&(z->deviations[zTerm]), &(z->deviations[zTerm]), &temp1, MPFR_RNDU);
        if (mpfr_cmp_ui(&(z->deviations[zTerm]), 0) < 0) {
            mpfr_set_ui(&(z->deviations[zTerm]), 0, MPFR_RNDN);
        }
    }
#endif // ARPRA_MIXED_TRIMMED_IAAA
    mpfi_intersect(&(z->true_range), &(z->true_range), &ia_range);
#endif // ARPRA_MIXED_IAAA

    // Handle domain violations.
    if (mpfr_nan_p(&(z->centre)) || mpfr_nan_p(&(z->radius))) {
        arpra_set_nan(z);
    }
    else if (mpfr_inf_p(&(z->centre)) || mpfr_inf_p(&(z->radius))) {
        arpra_set_inf(z);
    }

    // Clear vars.
    mpfr_clear(&temp1);
    mpfr_clear(&temp2);
    mpfi_clear(&ia_range);
}
