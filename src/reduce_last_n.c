/*
 * reduce_last_n.c -- Reduce the last n terms.
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

void mpfa_reduce_last_n (mpfa_ptr z, mpfa_uint_t n)
{
    mpfa_uint_t zTerm, zNext;
    mpfr_ptr *summands;
    mpfr_t temp;
    mpfa_prec_t prec_internal;

    // Handle trivial cases.
    if (n > z->nTerms) n = z->nTerms;
    if (n < 2) return;

    // Handle domain violations.
    if (mpfa_nan_p(z)) return;
    if (mpfa_inf_p(z)) return;

    // Initialise vars.
    prec_internal = mpfa_get_internal_prec();
    mpfr_init2(temp, prec_internal);
    mpfr_set_prec(&(z->radius), prec_internal);
    mpfr_set_ui(&(z->radius), 0, MPFR_RNDU);
    zTerm = z->nTerms - n;
    summands = malloc(n * sizeof(mpfr_ptr));

    // Add leading deviation terms to radius.
    for (zNext = 0; zNext < zTerm; zNext++) {
        mpfr_abs(temp, &(z->deviations[zNext]), MPFR_RNDU);
        mpfr_add(&(z->radius), &(z->radius), temp, MPFR_RNDU);
    }

    // Merge the last n deviation terms.
    for (zNext = zTerm; zNext < z->nTerms; zNext++) {
        mpfr_abs(&(z->deviations[zNext]), &(z->deviations[zNext]), MPFR_RNDN);
        summands[zNext - zTerm] = &(z->deviations[zNext]);
    }
    mpfr_sum(&(z->deviations[zTerm]), summands, n, MPFR_RNDU);

    // Store nonzero merged deviation term.
    if (!mpfr_zero_p(&(z->deviations[zTerm]))) {
        z->symbols[zTerm] = mpfa_next_sym();
        mpfr_add(&(z->radius), &(z->radius), &(z->deviations[zTerm]), MPFR_RNDU);
        zTerm++;
    }

    // Clear unused deviation terms.
    for (zNext = zTerm; zNext < z->nTerms; zNext++) {
        mpfr_clear(&(z->deviations[zNext]));
    }

    // Handle domain violations, and resize memory.
    z->nTerms = zTerm;
    if (mpfr_nan_p(&(z->centre)) || mpfr_nan_p(&(z->radius))) {
        mpfa_set_nan(z);
    }
    else if (mpfr_inf_p(&(z->centre)) || mpfr_inf_p(&(z->radius))) {
        mpfa_set_inf(z);
    }
    else {
        if (z->nTerms == 0) {
            free(z->symbols);
            free(z->deviations);
        }
        else {
            z->symbols = realloc(z->symbols, z->nTerms * sizeof(mpfa_uint_t));
            z->deviations = realloc(z->deviations, z->nTerms * sizeof(mpfr_t));
        }
    }

    // Clear vars.
    mpfr_clear(temp);
    free(summands);
}
