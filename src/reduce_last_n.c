/*
 * reduce_last_n.c -- Reduce the last n deviation terms.
 *
 * Copyright 2017-2018 James Paul Turner.
 *
 * This file is part of the Arpra library.
 *
 * The Arpra library is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published
 * by the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * The Arpra library is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
 * or FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public
 * License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with the Arpra library. If not, see <http://www.gnu.org/licenses/>.
 */

#include "arpra-impl.h"

void arpra_reduce_last_n (arpra_range *z, arpra_uint n)
{
    arpra_uint zTerm, zNext;
    arpra_mpfr **summands;
    arpra_mpfr z_lo, z_hi, temp;
    arpra_prec prec_internal;

    // Handle trivial cases.
    if (n > z->nTerms) n = z->nTerms;
    if (n < 2) return;

    // Handle domain violations.
    if (arpra_nan_p(z)) return;
    if (arpra_inf_p(z)) return;

    // Initialise vars.
    prec_internal = arpra_get_internal_precision();
    mpfr_init2(&z_lo, prec_internal);
    mpfr_init2(&z_hi, prec_internal);
    mpfr_init2(&temp, prec_internal);
    mpfr_set_ui(&(z->radius), 0, MPFR_RNDU);
    zTerm = z->nTerms - n;
    summands = malloc(n * sizeof(arpra_mpfr *));

    // Merge the last n deviation terms.
    for (zNext = zTerm; zNext < z->nTerms; zNext++) {
        mpfr_abs(&(z->deviations[zNext]), &(z->deviations[zNext]), MPFR_RNDN);
        summands[zNext - zTerm] = &(z->deviations[zNext]);
    }
    mpfr_sum(&(z->deviations[zTerm]), summands, n, MPFR_RNDU);

    // Add the remaining deviation terms to radius.
    for (zNext = 0; zNext < zTerm; zNext++) {
        mpfr_abs(&temp, &(z->deviations[zNext]), MPFR_RNDU);
        mpfr_add(&(z->radius), &(z->radius), &temp, MPFR_RNDU);
    }

    // Round the result to the target precision.
    mpfr_sub(&z_lo, &(z->centre), &(z->radius), MPFR_RNDD);
    mpfr_sub(&z_lo, &z_lo, &(z->deviations[zTerm]), MPFR_RNDD);
    if (mpfr_prec_round(&z_lo, z->precision, MPFR_RNDD)) {
        arpra_helper_error_ulp(&z_lo, &z_lo);
    }
    else {
        mpfr_set_ui(&z_lo, 0, MPFR_RNDN);
    }
    mpfr_add(&z_hi, &(z->centre), &(z->radius), MPFR_RNDU);
    mpfr_add(&z_hi, &z_hi, &(z->deviations[zTerm]), MPFR_RNDU);
    if (mpfr_prec_round(&z_hi, z->precision, MPFR_RNDU)) {
        arpra_helper_error_ulp(&z_hi, &z_hi);
    }
    else {
        mpfr_set_ui(&z_hi, 0, MPFR_RNDN);
    }
    mpfr_max(&temp, &z_lo, &z_hi, MPFR_RNDU);
    mpfr_add(&(z->deviations[zTerm]), &(z->deviations[zTerm]), &temp, MPFR_RNDU);

    // Store nonzero merged deviation term.
    if (!mpfr_zero_p(&(z->deviations[zTerm]))) {
        z->symbols[zTerm] = arpra_next_symbol();
        mpfr_add(&(z->radius), &(z->radius), &(z->deviations[zTerm]), MPFR_RNDU);
        zTerm++;
    }

    // Clear unused deviation terms.
    for (zNext = zTerm; zNext < z->nTerms; zNext++) {
        mpfr_clear(&(z->deviations[zNext]));
    }

    // Handle domain violations, and resize memory.
    z->nTerms = zTerm;
    if (mpfr_nan_p(&(z->radius))) {
        arpra_set_nan(z);
    }
    else if (mpfr_inf_p(&(z->radius))) {
        arpra_set_inf(z);
    }
    else {
        if (z->nTerms == 0) {
            free(z->symbols);
            free(z->deviations);
        }
        else {
            z->symbols = realloc(z->symbols, z->nTerms * sizeof(arpra_uint));
            z->deviations = realloc(z->deviations, z->nTerms * sizeof(arpra_mpfr));
        }
    }

    // Clear vars.
    mpfr_clear(&z_lo);
    mpfr_clear(&z_hi);
    mpfr_clear(&temp);
    free(summands);
}
