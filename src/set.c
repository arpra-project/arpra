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
    arpra_mpfr z_lo, z_hi, temp, error;
    arpra_prec prec_internal;

    // Handle trivial cases.
    if (z == x) return;

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
    mpfr_init2(&z_lo, prec_internal);
    mpfr_init2(&z_hi, prec_internal);
    mpfr_init2(&temp, prec_internal);
    mpfr_init2(&error, prec_internal);
    mpfr_set_prec(&(z->centre), prec_internal);
    mpfr_set_prec(&(z->radius), prec_internal);
    mpfr_set_ui(&error, 0, MPFR_RNDU);
    mpfr_set_ui(&(z->radius), 0, MPFR_RNDU);

    // z_0 = x_0
    if (mpfr_set(&(z->centre), &(x->centre), MPFR_RNDN)) {
        arpra_helper_error_half_ulp(&temp, &(z->centre));
        mpfr_add(&error, &error, &temp, MPFR_RNDU);
    }

    // Replace existing deviation term memory.
    arpra_clear_terms(z);
    z->nTerms = x->nTerms + 1;
    z->symbols = malloc(z->nTerms * sizeof(arpra_uint));
    z->deviations = malloc(z->nTerms * sizeof(arpra_range));

    // Copy deviation terms over.
    for (xTerm = 0, zTerm = 0; xTerm < x->nTerms; xTerm++) {
        z->symbols[zTerm] = x->symbols[zTerm];
        mpfr_init2(&(z->deviations[zTerm]), prec_internal);

        // z_i = x_i
        if (mpfr_set(&(z->deviations[zTerm]), &(x->deviations[xTerm]), MPFR_RNDN)) {
            arpra_helper_error_half_ulp(&temp, &(z->deviations[zTerm]));
            mpfr_add(&error, &error, &temp, MPFR_RNDU);
        }

        // Store nonzero deviation terms.
        if (mpfr_zero_p(&(z->deviations[zTerm]))) {
            mpfr_clear(&(z->deviations[zTerm]));
        }
        else {
            mpfr_abs(&temp, &(z->deviations[zTerm]), MPFR_RNDU);
            mpfr_add(&(z->radius), &(z->radius), &temp, MPFR_RNDU);
            zTerm++;
        }
    }

    // Round the result to the target precision.
    mpfr_sub(&z_lo, &(z->centre), &(z->radius), MPFR_RNDD);
    mpfr_sub(&z_lo, &z_lo, &error, MPFR_RNDD);
    if (mpfr_prec_round(&z_lo, z->precision, MPFR_RNDD)) {
        arpra_helper_error_ulp(&z_lo, &z_lo);
    }
    else {
        mpfr_set_ui(&z_lo, 0, MPFR_RNDN);
    }
    mpfr_add(&z_hi, &(z->centre), &(z->radius), MPFR_RNDU);
    mpfr_add(&z_hi, &z_hi, &error, MPFR_RNDU);
    if (mpfr_prec_round(&z_hi, z->precision, MPFR_RNDU)) {
        arpra_helper_error_ulp(&z_hi, &z_hi);
    }
    else {
        mpfr_set_ui(&z_hi, 0, MPFR_RNDN);
    }
    mpfr_max(&temp, &z_lo, &z_hi, MPFR_RNDU);
    mpfr_add(&error, &error, &temp, MPFR_RNDU);

    // Store nonzero numerical error term.
    if (!mpfr_zero_p(&error)) {
        z->symbols[zTerm] = arpra_next_symbol();
        mpfr_init2(&(z->deviations[zTerm]), prec_internal);
        mpfr_set(&(z->deviations[zTerm]), &error, MPFR_RNDU);
        mpfr_add(&(z->radius), &(z->radius), &(z->deviations[zTerm]), MPFR_RNDU);
        zTerm++;
    }

    // Handle domain violations, and resize memory.
    z->nTerms = zTerm;
    if (mpfr_nan_p(&(z->centre)) || mpfr_nan_p(&(z->radius))) {
        arpra_set_nan(z);
    }
    else if (mpfr_inf_p(&(z->centre)) || mpfr_inf_p(&(z->radius))) {
        arpra_set_inf(z);
    }
    else {
        if (z->nTerms == 0) {
            free(z->symbols);
            free(z->deviations);
        }
        else {
            z->symbols = realloc(z->symbols, z->nTerms * sizeof(arpra_uint));
            z->deviations = realloc(z->deviations, z->nTerms * sizeof(arpra_mpfr));
        }
    }

    // Clear vars.
    mpfr_clear(&z_lo);
    mpfr_clear(&z_hi);
    mpfr_clear(&temp);
    mpfr_clear(&error);
}
