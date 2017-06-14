/*
 * set_mpfr_rad.c -- Set centre and radius using MPFR numbers.
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

#include "mpfa.h"
#include <stdlib.h>

void mpfa_set_mpfr_rad (mpfa_ptr z, mpfr_srcptr centre, mpfr_srcptr radius) {
    unsigned zTerm;

    if (mpfr_set(&(z->centre), centre, MPFR_RNDN)) {
        mpfa_error(&(z->radius), &(z->centre));
        mpfr_add(&(z->radius), &(z->radius), radius, MPFR_RNDU);
    }
    else {
        mpfr_set(&(z->radius), radius, MPFR_RNDU);
    }

    if (mpfr_zero_p(&(z->radius))) {
        if (z->nTerms > 0) {
            for (zTerm = 0; zTerm < z->nTerms; zTerm++) {
                mpfr_clear(&(z->deviations[zTerm]));
            }
            z->nTerms = 0;
            free(z->symbols);
            free(z->deviations);
        }
    }
    else {
        if (z->nTerms == 0) {
            z->symbols = malloc(sizeof(unsigned));
            z->deviations = malloc(sizeof(mpfa_t));
            mpfr_init2(&(z->deviations[0]), mpfr_get_prec(&(z->centre)));
        }
        else if (z->nTerms >= 2) {
            for (zTerm = 1; zTerm < z->nTerms; zTerm++) {
                mpfr_clear(&(z->deviations[zTerm]));
            }
            z->symbols = realloc(z->symbols, sizeof(unsigned));
            z->deviations = realloc(z->deviations, sizeof(mpfa_t));
        }
        z->nTerms = 1;
        z->symbols[0] = mpfa_next_sym();
        mpfr_set(&(z->deviations[0]), &(z->radius), MPFR_RNDN);
    }
}
