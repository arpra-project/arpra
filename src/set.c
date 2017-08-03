/*
 * set.c -- Set one affine form with the values of another.
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

void mpfa_set (mpfa_ptr z, mpfa_srcptr x) {
    mpfa_uint_t xTerm, zTerm;
    mpfa_prec_t prec;
    mpfr_t temp, error;

    if (z == x) return;

    prec = mpfr_get_prec(&(z->centre));
    mpfr_inits2(prec, temp, error, (mpfr_ptr) NULL);
    mpfr_set_si(error, 0, MPFR_RNDN);
    mpfr_set_si(&(z->radius), 0, MPFR_RNDN);

    if (mpfr_set(&(z->centre), &(x->centre), MPFR_RNDN)) {
        mpfa_error(temp, &(z->centre));
        mpfr_add(error, error, temp, MPFR_RNDU);
    }

    if (z->nTerms < x->nTerms) {
        // need to grow z and initialise new terms
        if (z->nTerms == 0) {
            z->symbols = malloc((x->nTerms + 1) * sizeof(mpfa_uint_t));
            z->deviations = malloc((x->nTerms + 1) * sizeof(mpfa_t));
        }
        else {
            z->symbols = realloc(z->symbols, (x->nTerms + 1) * sizeof(mpfa_uint_t));
            z->deviations = realloc(z->deviations, (x->nTerms + 1) * sizeof(mpfa_t));
        }
        for (zTerm = z->nTerms; zTerm < x->nTerms; zTerm++) {
            mpfr_init2(&(z->deviations[zTerm]), prec);
        }
    }
    else if (z->nTerms > x->nTerms) {
        // need to clear unused terms
        for (zTerm = x->nTerms; zTerm < z->nTerms; zTerm++) {
            mpfr_clear(&(z->deviations[zTerm]));
        }
    }

    for (xTerm = 0, zTerm = 0; xTerm < x->nTerms; xTerm++) {
        z->symbols[zTerm] = x->symbols[zTerm];

        if (mpfr_set(&(z->deviations[zTerm]), &(x->deviations[xTerm]), MPFR_RNDN)) {
            mpfa_error(temp, &(z->deviations[zTerm]));
            mpfr_add(error, error, temp, MPFR_RNDU);
        }

        if (mpfr_zero_p(&(z->deviations[zTerm]))) {
            mpfr_clear(&(z->deviations[zTerm]));
        }
        else {
            mpfr_abs(temp, &(z->deviations[zTerm]), MPFR_RNDN);
            mpfr_add(&(z->radius), &(z->radius), temp, MPFR_RNDU);
            zTerm++;
        }
    }

    if (!mpfr_zero_p(error)) {
        z->symbols[zTerm] = mpfa_next_sym();
        mpfr_init2(&(z->deviations[zTerm]), prec);
        mpfr_set(&(z->deviations[zTerm]), error, MPFR_RNDN);
        mpfr_add(&(z->radius), &(z->radius), error, MPFR_RNDU);
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

    mpfr_clears(temp, error, (mpfr_ptr) NULL);
}
