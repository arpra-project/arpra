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

#include "mpfa-impl.h"

void mpfa_condense_small (mpfa_ptr z, double fraction)
{
    mpfa_uint_t zTerm, zNext;
    mpfr_t temp, error, threshold;
    mpfa_prec_t prec, prec_internal;

    // Check input, and handle trivial case.
    if ((z->nTerms < 2) || (fraction >= 1)) return;

    // Init temp vars, and set internal precision.
    prec = mpfa_get_prec(z);
    prec_internal = mpfa_get_internal_prec();
    mpfr_init2(temp, prec_internal);
    mpfr_init2(error, prec_internal);
    mpfr_init2(threshold, prec_internal);
    mpfr_set_si(error, 0, MPFR_RNDN);
    mpfr_mul_d(threshold, &(z->radius), fraction, MPFR_RNDN);
    mpfr_prec_round(&(z->radius), prec_internal, MPFR_RNDU);
    mpfr_set_si(&(z->radius), 0, MPFR_RNDN);
    zTerm = 0;

    for (zNext = 0; zNext < z->nTerms; zNext++) {
        mpfr_abs(temp, &(z->deviations[zNext]), MPFR_RNDN);

        if (mpfr_lessequal_p(temp, threshold)) {
            // If noise term is smaller than threshold, condense it.
            mpfr_add(error, error, temp, MPFR_RNDU);
        }
        else {
            // Else shift noise term up, and add it to radius.
            mpfr_prec_round(&(z->deviations[zTerm]), prec_internal, MPFR_RNDN);
            if (zTerm < zNext) {
                z->symbols[zTerm] = z->symbols[zNext];
                mpfr_set(&(z->deviations[zTerm]), &(z->deviations[zNext]), MPFR_RNDN);
            }
            mpfr_add(&(z->radius), &(z->radius), temp, MPFR_RNDU);
            zTerm++;
        }
    }

    // Store nonzero condensed noise term, and add it to radius.
    if (!mpfr_zero_p(error)) {
        z->symbols[zTerm] = mpfa_next_sym();
        mpfr_prec_round(&(z->deviations[zTerm]), prec_internal, MPFR_RNDN);
        mpfr_set(&(z->deviations[zTerm]), error, MPFR_RNDU);
        mpfr_add(&(z->radius), &(z->radius), &(z->deviations[zTerm]), MPFR_RNDU);
        zTerm++;
    }

    // Clear unused noise terms.
    for (zNext = zTerm; zNext < z->nTerms; zNext++) {
        mpfr_clear(&(z->deviations[zNext]));
    }

    // Resize noise term memory.
    if (zTerm == 0) {
        free(z->symbols);
        free(z->deviations);
    }
    else {
        z->symbols = realloc(z->symbols, zTerm * sizeof(mpfa_uint_t));
        z->deviations = realloc(z->deviations, zTerm * sizeof(mpfr_t));
    }
    z->nTerms = zTerm;

    // Clear temp vars.
    mpfr_clear(temp);
    mpfr_clear(error);
    mpfr_clear(threshold);
}
