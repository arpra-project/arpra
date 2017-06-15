/*
 * condense_last_n.c -- Condense the last n deviation terms into one.
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

#include "mpfa.h"
#include <stdlib.h>

void mpfa_condense_last_n (mpfa_ptr z, mpfa_uint_t n) {
    mpfa_uint_t zTerm, last;

    if (n < 2) return;

    mpfr_set_si(&(z->radius), 0, MPFR_RNDN);
    last = z->nTerms - n;

    for (zTerm = (last + 1); zTerm < z->nTerms; zTerm++) {
        mpfr_add(&(z->deviations[last]), &(z->deviations[last]), &(z->deviations[zTerm]), MPFR_RNDU);
        mpfr_clear(&(z->deviations[zTerm]));
    }

    z->nTerms = last + 1;
    z->symbols = realloc(z->symbols, z->nTerms * sizeof(mpfa_uint_t));
    z->deviations = realloc(z->deviations, z->nTerms * sizeof(mpfa_t));

    for (zTerm = 0; zTerm < z->nTerms; zTerm++) {
        mpfr_add(&(z->radius), &(z->radius), &(z->deviations[zTerm]), MPFR_RNDU);
    }
}
