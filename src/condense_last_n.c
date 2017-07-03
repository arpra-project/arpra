/*
 * condense_last_n.c -- Condense the last n terms.
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
    mpfa_uint_t zTerm, zNext;
    mpfr_prec_t prec;
    mpfr_t temp;

    if ((z->nTerms < 2) || (n < 2)) return;

    if (n > z->nTerms) {
        zTerm = 0
    }
    else {
        zTerm = z->nTerms - n;
    }

    prec = mpfr_get_prec(&(z->centre));
    mpfr_init2(temp, prec);
    mpfr_set_si(&(z->radius), 0, MPFR_RNDN);

    for (zNext = (zTerm + 1); zNext < z->nTerms; zNext++) {
        mpfr_add(&(z->deviations[zTerm]), &(z->deviations[zTerm]), &(z->deviations[zNext]), MPFR_RNDU);
        mpfr_clear(&(z->deviations[zNext]));
    }

    if (mpfr_zero_p(&(z->deviations[zTerm]))) {
        mpfr_clear(&(z->deviations[zTerm]));
    }
    else {
        z->symbols[zTerm] = mpfa_next_sym();
        zTerm++;
    }

    if (z->nTerms > 0) {
        if (zTerm == 0) {
            free(z->symbols);
            free(z->deviations);
        }
        else {
            z->symbols = realloc(z->symbols, zTerm * sizeof(mpfa_uint_t));
            z->deviations = realloc(z->deviations, zTerm * sizeof(mpfr_t));
        }
    }
    z->nTerms = zTerm;

    for (zTerm = 0; zTerm < z->nTerms; zTerm++) {
        mpfr_abs(temp, &(z->deviations[zTerm]), MPFR_RNDN);
        mpfr_add(&(z->radius), &(z->radius), temp, MPFR_RNDU);
    }

    mpfr_clear(temp);
}
