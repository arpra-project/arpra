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

void mpfa_set_str (mpfa_ptr z, const char *centre, mpfa_int_t base) {
    mpfa_uint_t zTerm;
    mpfa_prec_t prec_internal;

    prec_internal = mpfa_get_internal_prec();

    // If centre has rounding error:
    if (mpfr_set_str(&(z->centre), centre, base, MPFR_RNDN)) {
        mpfa_error(&(z->radius), &(z->centre));
    }

    // Else centre has no rounding error:
    else {
        mpfr_set_si(&(z->radius), 0, MPFR_RNDU);
    }

    // Clear existing noise terms.
    mpfa_clear_terms(z);

    // If radius is nonzero:
    if (!mpfr_zero_p(&(z->radius))) {
        // Allocate one noise term.
        z->nTerms = 1;
        z->symbols = malloc(sizeof(mpfa_uint_t));
        z->deviations = malloc(sizeof(mpfa_t));

        // Set noise term.
        z->symbols[0] = mpfa_next_sym();
        mpfr_init2(&(z->deviations[0]), prec_internal);
        mpfr_set(&(z->deviations[0]), &(z->radius), MPFR_RNDN);
    }
}
