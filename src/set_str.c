/*
 * set_str.c -- Set centre using a C string.
 *
 * Copyright 2017 James Paul Turner.
 *
 * This file is part of the MPFA library.
 *
 * The MPFA library is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published
 * by the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * The MPFA library is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
 * or FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public
 * License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with the MPFA library. If not, see <http://www.gnu.org/licenses/>.
 */

#include "mpfa-impl.h"

void mpfa_set_str (mpfa_ptr z, const char *centre, mpfa_int_t base)
{
    mpfa_prec_t prec, prec_internal;

    // Initialise vars.
    prec = mpfa_get_prec(z);
    prec_internal = mpfa_get_internal_prec();
    mpfr_prec_round(&(z->radius), prec_internal, MPFR_RNDU);
    mpfr_set_ui(&(z->radius), 0, MPFR_RNDU);

    // Add centre rounding error to deviation.
    if (mpfr_set_str(&(z->centre), centre, base, MPFR_RNDN)) {
        mpfa_error(&(z->radius), &(z->centre));
    }

    // Clear existing deviation terms.
    mpfa_clear_terms(z);

    // If radius is nonzero:
    if (!mpfr_zero_p(&(z->radius))) {
        // Allocate one deviation term.
        z->nTerms = 1;
        z->symbols = malloc(sizeof(mpfa_uint_t));
        z->deviations = malloc(sizeof(mpfa_t));

        // Set deviation term.
        z->symbols[0] = mpfa_next_sym();
        mpfr_init2(&(z->deviations[0]), prec);
        mpfr_set(&(z->deviations[0]), &(z->radius), MPFR_RNDU);
        mpfr_set(&(z->radius), &(z->deviations[0]), MPFR_RNDU);
    }

    // Handle domain violations.
    if (mpfr_nan_p(&(z->centre)) || mpfr_nan_p(&(z->radius))) {
        mpfa_set_nan(z);
    }
    else if (mpfr_inf_p(&(z->centre)) || mpfr_inf_p(&(z->radius))) {
        mpfa_set_inf(z);
    }
}
