/*
 * condense_small.c -- Condense terms less than some fraction of the radius.
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

void mpfa_condense_small (mpfa_ptr z, double fraction) {
    mpfa_uint_t zTerm, zNext;
    mpfr_prec_t prec;
    mpfr_t temp, error, threshold;

    if (z->nTerms < 2) return;

    // check double fraction in [0, 1]






    prec = mpfr_get_prec(&(z->centre));
    mpfr_inits2(prec, temp, error, threshold, (mpfr_ptr) NULL);
    mpfr_mul_d(threshold, &(z->radius), fraction, MPFR_RNDN);
    mpfr_set_si(error, 0, MPFR_RNDN);
    mpfr_set_si(&(z->radius), 0, MPFR_RNDN);

    for (zTerm = 0, zNext = 0; zNext < z->nTerms; zNext++) {
        mpfr_abs(temp, &(z->deviations[zTerm]), MPFR_RNDN);

        if (mpfr_lessequal_p(temp, threshold)) {
            
        }
        else {
            
        }
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



    mpfr_clears(temp, error, threshold, (mpfr_ptr) NULL);
}
