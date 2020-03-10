/*
 * set_str_rad.c -- Set centre and radius using C strings.
 *
 * Copyright 2017-2018 James Paul Turner.
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

void arpra_set_str_rad (arpra_range *z, const char *centre, const char *radius, const arpra_int base)
{
    arpra_mpfr temp;
    arpra_prec prec_internal;

    // Initialise vars.
    prec_internal = arpra_get_internal_precision();
    mpfr_init2(&temp, prec_internal);
    mpfr_set_prec(&(z->centre), prec_internal);
    mpfr_set_prec(&(z->radius), prec_internal);

    // Handle negative radius.
    mpfr_set_str(&temp, radius, base, MPFR_RNDU);
    if (mpfr_cmp_ui(&temp, 0) < 0) {
        arpra_set_nan(z);
        mpfr_clear(&temp);
        return;
    }

    // MPFI set_str
    mpfr_set_str(&(z->true_range.left), centre, base, MPFR_RNDD);
    mpfr_set_str(&(z->true_range.right), centre, base, MPFR_RNDU);
    mpfr_sub(&(z->true_range.left), &(z->true_range.left), &temp, MPFR_RNDD);
    mpfr_add(&(z->true_range.right), &(z->true_range.right), &temp, MPFR_RNDU);

    // z_0 = (x_lo + x_hi) / 2
    mpfi_mid(&(z->centre), &(z->true_range));

    // rad(z) = max{(z_0 - x_lo), (x_hi - z_0)}
    mpfr_sub(&(z->radius), &(z->centre), &(z->true_range.left), MPFR_RNDU);
    mpfr_sub(&temp, &(z->true_range.right), &(z->centre), MPFR_RNDU);
    mpfr_max(&(z->radius), &(z->radius), &temp, MPFR_RNDU);

    // Clear existing deviation terms.
    arpra_helper_clear_terms(z);

    // Store nonzero rounding error term.
    if (!mpfr_zero_p(&(z->radius))) {
        z->nTerms = 1;
        z->symbols = malloc(sizeof(arpra_uint));
        z->deviations = malloc(sizeof(arpra_mpfr));
        z->symbols[0] = arpra_helper_next_symbol();
        mpfr_init2(&(z->deviations[0]), prec_internal);
        mpfr_set(&(z->deviations[0]), &(z->radius), MPFR_RNDU);
    }

    // Check for NaN and Inf.
    arpra_helper_check_result(z);

    // Clear vars.
    mpfr_clear(&temp);
}
