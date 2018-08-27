/*
 * reduce_small.c -- Reduce deviation terms less than a given threshold.
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

void arpra_reduce_small (arpra_range *z, double min_fraction)
{
    arpra_uint zTerm, zNext;
    arpra_mpfr threshold;
    arpra_precision prec_internal;

    // Handle trivial cases.
    if (z->nTerms < 2) return;
    if (min_fraction <= 0.0) return;

    // Handle domain violations.
    if (!arpra_bounded_p(z)) return;

    // Initialise vars.
    prec_internal = arpra_get_internal_precision();
    mpfr_init2(&threshold, prec_internal);
    mpfr_mul_d(&threshold, &(z->radius), min_fraction, MPFR_RNDN);

    // Move small deviation terms to back.
    for (zNext = 0, zTerm = 0; zNext < z->nTerms; zNext++) {
        if (mpfr_cmpabs(&(z->deviations[zNext]), &threshold) >= 0) {
            if (zTerm < zNext) {
                z->symbols[zTerm] = z->symbols[zNext];
                mpfr_swap(&(z->deviations[zTerm]), &(z->deviations[zNext]));
            }
            zTerm++;
        }
    }

    // Merge the small deviation terms.
    arpra_reduce_last_n(z, (z->nTerms - zTerm));

    // Clear vars.
    mpfr_clear(&threshold);
}

void arpra_reduce_small_abs (arpra_range *z, const arpra_mpfr *abs_threshold)
{
    arpra_uint zTerm, zNext;

    // Handle trivial cases.
    if (z->nTerms < 2) return;
    if (mpfr_cmp_ui(abs_threshold, 0) <= 0) return;

    // Handle domain violations.
    if (!arpra_bounded_p(z)) return;

    // Move small deviation terms to the back.
    for (zNext = 0, zTerm = 0; zNext < z->nTerms; zNext++) {
        if (mpfr_cmpabs(&(z->deviations[zNext]), abs_threshold) >= 0) {
            if (zTerm < zNext) {
                z->symbols[zTerm] = z->symbols[zNext];
                mpfr_swap(&(z->deviations[zTerm]), &(z->deviations[zNext]));
            }
            zTerm++;
        }
    }

    // Reduce small deviation terms.
    arpra_reduce_last_n(z, (z->nTerms - zTerm));
}

void arpra_reduce_small_rel (arpra_range *z, const arpra_mpfr *rel_threshold)
{
    arpra_mpfr abs_threshold;
    arpra_precision prec_internal;

    // Handle trivial cases.
    if (z->nTerms < 2) return;
    if (mpfr_cmp_ui(rel_threshold, 0) <= 0) return;

    // Handle domain violations.
    if (!arpra_bounded_p(z)) return;

    // Initialise vars.
    prec_internal = arpra_get_internal_precision();
    mpfr_init2(&abs_threshold, prec_internal);

    // Convert to absolute threshold.
    mpfr_mul(&abs_threshold, &(z->radius), rel_threshold, MPFR_RNDU);
    arpra_reduce_small_abs(z, &abs_threshold);

    // Clear vars.
    mpfr_clear(&abs_threshold);
}
