/*
 * set_mpfr.c -- Set centre using an MPFR number.
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

void arpra_set_mpfr (arpra_range *z, const arpra_mpfr *centre)
{
    arpra_mpfr temp;
    arpra_mpfi range;
    arpra_prec prec_internal;

    // Handle domain violations.
    if (mpfr_nan_p(centre)) {
        arpra_set_nan(z);
        return;
    }
    if (mpfr_inf_p(centre)) {
        arpra_set_inf(z);
        return;
    }

    // Initialise vars.
    prec_internal = arpra_get_internal_precision();
    mpfr_init2(&temp, prec_internal);
    mpfi_init2(&range, prec_internal);
    mpfr_set_prec(&(z->centre), prec_internal);
    mpfr_set_prec(&(z->radius), prec_internal);
    mpfr_set_ui(&(z->radius), 0, MPFR_RNDU);

    // Add centre rounding error to deviation.
    if (mpfr_set(&(z->centre), centre, MPFR_RNDN)) {
        arpra_helper_error_half_ulp(&temp, &(z->centre));
        mpfr_add(&(z->radius), &(z->radius), &temp, MPFR_RNDU);
    }

    // Clear existing deviation terms.
    arpra_clear_terms(z);

    // Round the result to the target precision.
    mpfr_sub(&(range.left), &(z->centre), &(z->radius), MPFR_RNDD);
    if (mpfr_prec_round(&(range.left), z->precision, MPFR_RNDD)) {
        arpra_helper_error_ulp(&(range.left), &(range.left));
    }
    else {
        mpfr_set_ui(&(range.left), 0, MPFR_RNDN);
    }
    mpfr_add(&(range.right), &(z->centre), &(z->radius), MPFR_RNDU);
    if (mpfr_prec_round(&(range.right), z->precision, MPFR_RNDU)) {
        arpra_helper_error_ulp(&(range.right), &(range.right));
    }
    else {
        mpfr_set_ui(&(range.right), 0, MPFR_RNDN);
    }
    mpfr_max(&temp, &(range.left), &(range.right), MPFR_RNDU);
    mpfr_add(&(z->radius), &(z->radius), &temp, MPFR_RNDU);

    // Store nonzero numerical error term.
    if (!mpfr_zero_p(&(z->radius))) {
        z->nTerms = 1;
        z->symbols = malloc(sizeof(arpra_uint));
        z->deviations = malloc(sizeof(arpra_range));
        z->symbols[0] = arpra_next_symbol();
        mpfr_init2(&(z->deviations[0]), prec_internal);
        mpfr_set(&(z->deviations[0]), &(z->radius), MPFR_RNDU);
    }

    // Compute true range.
    mpfr_sub(&(z->true_range.left), &(z->centre), &(z->radius), MPFR_RNDD);
    mpfr_add(&(z->true_range.right), &(z->centre), &(z->radius), MPFR_RNDU);

    // Handle domain violations.
    if (mpfr_nan_p(&(z->centre)) || mpfr_nan_p(&(z->radius))) {
        arpra_set_nan(z);
    }
    else if (mpfr_inf_p(&(z->centre)) || mpfr_inf_p(&(z->radius))) {
        arpra_set_inf(z);
    }

    // Clear vars.
    mpfr_clear(&temp);
    mpfi_clear(&range);
}
