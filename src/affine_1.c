/*
 * affine_1.c -- Compute a univariate affine function of an arpra_t.
 *
 * Copyright 2016-2017 James Paul Turner.
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

void arpra_affine_1 (arpra_ptr z, arpra_srcptr x, mpfr_srcptr alpha, mpfr_srcptr gamma, mpfr_srcptr delta)
{
    arpra_uint_t xTerm, zTerm;
    mpfr_t temp, error;
    arpra_prec_t prec, prec_internal;
    arpra_t zNew;

    // Domain violations:
    // NaN  =  Nan
    // Inf  =  Inf

    // Handle domain violations.
    if (arpra_nan_p(x)) {
        arpra_set_nan(z);
        return;
    }
    if (arpra_inf_p(x)) {
        arpra_set_inf(z);
        return;
    }

    // Initialise vars.
    prec = arpra_get_prec(z);
    prec_internal = arpra_get_internal_prec();
    mpfr_init2(temp, prec_internal);
    mpfr_init2(error, prec_internal);
    mpfr_init2(&(zNew->centre), prec);
    mpfr_init2(&(zNew->radius), prec_internal);
    mpfr_set(error, delta, MPFR_RNDU);
    mpfr_set_si(&(zNew->radius), 0, MPFR_RNDU);

    // z_0 = (alpha * x_0) + gamma
    if (mpfr_fma(&(zNew->centre), alpha, &(x->centre), gamma, MPFR_RNDN)) {
        arpra_error(temp, &(zNew->centre));
        mpfr_add(error, error, temp, MPFR_RNDU);
    }

    // Allocate memory for all possible deviation terms.
    zNew->nTerms = x->nTerms + 1;
    zNew->symbols = malloc(zNew->nTerms * sizeof(arpra_uint_t));
    zNew->deviations = malloc(zNew->nTerms * sizeof(mpfr_t));

    for (xTerm = 0, zTerm = 0; xTerm < x->nTerms; xTerm++) {
        zNew->symbols[zTerm] = x->symbols[xTerm];
        mpfr_init2(&(zNew->deviations[zTerm]), prec);

        // z_i = (alpha * x_i)
        if (mpfr_mul(&(zNew->deviations[zTerm]), alpha, &(x->deviations[xTerm]), MPFR_RNDN)) {
            arpra_error(temp, &(zNew->deviations[zTerm]));
            mpfr_add(error, error, temp, MPFR_RNDU);
        }

        // Store nonzero deviation terms.
        if (mpfr_zero_p(&(zNew->deviations[zTerm]))) {
            mpfr_clear(&(zNew->deviations[zTerm]));
        }
        else {
            mpfr_abs(temp, &(zNew->deviations[zTerm]), MPFR_RNDU);
            mpfr_add(&(zNew->radius), &(zNew->radius), temp, MPFR_RNDU);
            zTerm++;
        }
    }

    // Store nonzero numerical error term.
    if (!mpfr_zero_p(error)) {
        zNew->symbols[zTerm] = arpra_next_sym();
        mpfr_init2(&(zNew->deviations[zTerm]), prec);
        mpfr_set(&(zNew->deviations[zTerm]), error, MPFR_RNDU);
        mpfr_add(&(zNew->radius), &(zNew->radius), &(zNew->deviations[zTerm]), MPFR_RNDU);
        zTerm++;
    }

    // Handle domain violations, and free unused memory.
    zNew->nTerms = zTerm;
    if (mpfr_nan_p(&(zNew->centre)) || mpfr_nan_p(&(zNew->radius))) {
        arpra_set_nan(zNew);
    }
    else if (mpfr_inf_p(&(zNew->centre)) || mpfr_inf_p(&(zNew->radius))) {
        arpra_set_inf(zNew);
    }
    else {
        if (zNew->nTerms == 0) {
            free(zNew->symbols);
            free(zNew->deviations);
        }
    }

    // Clear vars, and set z.
    mpfr_clear(temp);
    mpfr_clear(error);
    arpra_set(z, zNew);
    arpra_clear(zNew);
}
