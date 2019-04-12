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
    arpra_mpfr temp1, temp2;
    arpra_prec prec_internal;

    // Handle trivial cases.
    if ((z->nTerms == 0) || (n == 0)) return;
    if (n > z->nTerms) n = z->nTerms;

    // Handle domain violations.
    if (arpra_nan_p(z)) return;
    if (arpra_inf_p(z)) return;

    // Initialise vars.
    prec_internal = arpra_get_internal_precision();
    mpfr_init2(&temp1, prec_internal + 8);
    mpfr_init2(&temp2, prec_internal + 8);
    mpfr_set_ui(&(z->radius), 0, MPFR_RNDU);
    zTerm = z->nTerms - n;
    summands = malloc(n * sizeof(arpra_mpfr *));

    // Merge the last n deviation terms.
    for (zNext = zTerm; zNext < z->nTerms; zNext++) {
        mpfr_abs(&(z->deviations[zNext]), &(z->deviations[zNext]), MPFR_RNDN);
        summands[zNext - zTerm] = &(z->deviations[zNext]);
    }
    z->symbols[zTerm] = arpra_next_symbol();
    mpfr_sum(&(z->deviations[zTerm]), summands, n, MPFR_RNDU);

    // Clear the unused deviation terms.
    for (zNext = zTerm + 1; zNext < z->nTerms; zNext++) {
        mpfr_clear(&(z->deviations[zNext]));
    }

    // Add the remaining deviation terms to radius.
    for (zNext = 0; zNext < zTerm; zNext++) {
        mpfr_abs(&temp1, &(z->deviations[zNext]), MPFR_RNDU);
        mpfr_add(&(z->radius), &(z->radius), &temp1, MPFR_RNDU);
    }
    mpfr_add(&(z->radius), &(z->radius), &(z->deviations[zTerm]), MPFR_RNDU);
    z->nTerms = zTerm + 1;

#ifdef ARPRA_MIXED_TRIMMED_IAAA
    // Trim error term if AA range fully encloses mixed IA/AA range.
    mpfr_sub(&temp1, &(z->centre), &(z->radius), MPFR_RNDD);
    mpfr_add(&temp2, &(z->centre), &(z->radius), MPFR_RNDU);
    if (mpfr_less_p(&temp1, &(z->true_range.left))
        && mpfr_greater_p(&temp2, &(z->true_range.right))) {
        mpfr_sub(&temp1, &(z->true_range.left), &temp1, MPFR_RNDD);
        mpfr_sub(&temp2, &temp2, &(z->true_range.right), MPFR_RNDD);
        mpfr_min(&temp1, &temp1, &temp2, MPFR_RNDD);
        if (mpfr_greater_p(&temp1, &(z->deviations[z->nTerms - 1]))) {
            mpfr_sub(&(z->radius), &(z->radius), &(z->deviations[z->nTerms - 1]), MPFR_RNDU);
            mpfr_set_ui(&(z->deviations[z->nTerms - 1]), 0, MPFR_RNDZ);
        }
        else {
            mpfr_sub(&(z->radius), &(z->radius), &temp1, MPFR_RNDU);
            mpfr_sub(&(z->deviations[z->nTerms - 1]), &(z->deviations[z->nTerms - 1]), &temp1, MPFR_RNDU);
        }
    }
#endif // ARPRA_MIXED_TRIMMED_IAAA

    // Handle domain violations.
    if (mpfr_nan_p(&(z->radius))) {
        arpra_set_nan(z);
    }
    else if (mpfr_inf_p(&(z->radius))) {
        arpra_set_inf(z);
    }

    // Clear vars.
    mpfr_clear(&temp1);
    mpfr_clear(&temp2);
    free(summands);
}
