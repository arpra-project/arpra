/*
 * set_d_rad.c -- Set centre and radius using double-precision floats.
 *
 * Copyright 2016-2017 James Paul Turner.
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

void mpfa_set_d_rad (mpfa_ptr z, const double centre, const double radius) {
    mpfa_uint_t zTerm;
    mpfa_prec_t prec, prec_internal;

    prec = mpfa_get_prec(z);
    prec_internal = mpfa_get_internal_prec();

    // If centre has rounding error:
    if (mpfr_set_d(&(z->centre), centre, MPFR_RNDN)) {
        mpfr_prec_round(&(z->radius), prec_internal, MPFR_RNDU);
        mpfa_error(&(z->radius), &(z->centre));
        mpfr_add_d(&(z->radius), &(z->radius), radius, MPFR_RNDU);
        mpfr_prec_round(&(z->radius), prec, MPFR_RNDU);
    }

    // Else centre has no rounding error:
    else {
        mpfr_set_d(&(z->radius), radius, MPFR_RNDU);
    }

    // If radius is zero:
    if (mpfr_zero_p(&(z->radius))) {
        // Clear noise terms.
        if (z->nTerms > 0) {
            for (zTerm = 0; zTerm < z->nTerms; zTerm++) {
                mpfr_clear(&(z->deviations[zTerm]));
            }
            z->nTerms = 0;
            free(z->symbols);
            free(z->deviations);
        }
    }

    // Else radius is nonzero:
    else {
        if (z->nTerms == 0) {
            // If z has no noise term memory, allocate enough for one term.
            z->symbols = malloc(sizeof(mpfa_uint_t));
            z->deviations = malloc(sizeof(mpfa_t));
            mpfr_init2(&(z->deviations[0]), prec_internal);
        }
        else if (z->nTerms >= 2) {
            // Else if z has too much noise term memory, clear unused terms.
            mpfr_prec_round(&(z->deviations[0]), prec_internal, MPFR_RNDN);
            for (zTerm = 1; zTerm < z->nTerms; zTerm++) {
                mpfr_clear(&(z->deviations[zTerm]));
            }
            z->symbols = realloc(z->symbols, sizeof(mpfa_uint_t));
            z->deviations = realloc(z->deviations, sizeof(mpfa_t));
        }
        z->nTerms = 1;

        // Set noise term.
        z->symbols[0] = mpfa_next_sym();
        mpfr_set(&(z->deviations[0]), &(z->radius), MPFR_RNDN);
    }
}
