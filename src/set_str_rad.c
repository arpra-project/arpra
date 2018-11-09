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
    arpra_prec prec, prec_internal;
    arpra_mpfr temp;

    // Initialise vars.
    prec = arpra_get_precision(z);
    prec_internal = arpra_get_internal_precision();
    mpfr_init2(&temp, prec_internal);
    mpfr_set_str(&temp, radius, base, MPFR_RNDU);
    mpfr_set_prec(&(z->radius), prec_internal);
    mpfr_abs(&(z->radius), &temp, MPFR_RNDU);

    // Add centre rounding error to deviation.
    if (mpfr_set_str(&(z->centre), centre, base, MPFR_RNDN)) {
        arpra_helper_error_half_ulp(&temp, &(z->centre));
        mpfr_add(&(z->radius), &(z->radius), &temp, MPFR_RNDU);
    }

    // Clear existing deviaion terms.
    arpra_clear_terms(z);

    // If radius is nonzero:
    if (!mpfr_zero_p(&(z->radius))) {
        // Allocate one deviation term.
        z->nTerms = 1;
        z->symbols = malloc(sizeof(arpra_uint));
        z->deviations = malloc(sizeof(arpra_range));

        // Set deviation term.
        z->symbols[0] = arpra_next_symbol();
        mpfr_init2(&(z->deviations[0]), prec);
        mpfr_set(&(z->deviations[0]), &(z->radius), MPFR_RNDU);
        mpfr_set(&(z->radius), &(z->deviations[0]), MPFR_RNDU);
    }

    // Handle domain violations.
    if (mpfr_nan_p(&(z->centre)) || mpfr_nan_p(&(z->radius))) {
        arpra_set_nan(z);
    }
    else if (mpfr_inf_p(&(z->centre)) || mpfr_inf_p(&(z->radius))) {
        arpra_set_inf(z);
    }

    // Clear vars.
    mpfr_clear(&temp);
}
