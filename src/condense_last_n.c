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

#include "mpfa-impl.h"

void mpfa_condense_last_n (mpfa_ptr z, mpfa_uint_t n) {
    mpfa_uint_t zTerm, zNext;
    mpfr_ptr *summands;
    mpfr_t temp;
    mpfa_prec_t prec, prec_internal;

    // Check input, and handle trivial case.
    if (n > z->nTerms) n = z->nTerms;
    if (n < 2) return;

    // Init temp vars, and set internal precision.
    prec = mpfa_get_prec(z);
    prec_internal = mpfa_get_internal_prec();
    mpfr_init2(temp, prec_internal);
    mpfr_prec_round(&(z->radius), prec_internal, MPFR_RNDU);
    mpfr_set_si(&(z->radius), 0, MPFR_RNDN);
    zTerm = z->nTerms - n;
    summands = malloc(n * sizeof(mpfr_ptr));

    // Add leading noise terms to radius.
    for (zNext = 0; zNext < zTerm; zNext++) {
        mpfr_prec_round(&(z->deviations[zNext]), prec_internal, MPFR_RNDN);
        mpfr_abs(temp, &(z->deviations[zNext]), MPFR_RNDN);
        mpfr_add(&(z->radius), &(z->radius), temp, MPFR_RNDU);
    }

    // Condense the last n noise terms.
    for (zNext = zTerm; zNext < z->nTerms; zNext++) {
        mpfr_abs(&(z->deviations[zNext]), &(z->deviations[zNext]), MPFR_RNDN);
        summands[zNext - zTerm] = &(z->deviations[zNext]);
    }
    mpfr_prec_round(&(z->deviations[zTerm]), prec_internal, MPFR_RNDN);
    mpfr_sum(&(z->deviations[zTerm]), summands, n, MPFR_RNDU);

    // Store nonzero condensed noise term, and add it to radius.
    if (!mpfr_zero_p(&(z->deviations[zTerm]))) {
        z->symbols[zTerm] = mpfa_next_sym();
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

    // Round internal precision of radius to working precision.
    mpfr_prec_round(&(z->radius), prec, MPFR_RNDU);

    // Clear temp vars.
    mpfr_clear(temp);
    free(summands);
}
