/*
 * affine_1.c -- Compute a univariate affine function of an Arpra range.
 *
 * Copyright 2016-2018 James Paul Turner.
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

void arpra_affine_1 (arpra_range *z, const arpra_range *x,
                     const arpra_mpfr *alpha, const arpra_mpfr *gamma,
                     const arpra_mpfr *delta)
{
    arpra_uint xTerm, zTerm;
    arpra_mpfr temp1, temp2, error;
    arpra_range zNew;
    arpra_prec prec_internal;

    // Domain violations:
    // NaN  =  NaN
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
    prec_internal = arpra_get_internal_precision();
    mpfr_init2(&temp1, prec_internal + 8);
    mpfr_init2(&temp2, prec_internal + 8);
    mpfr_init2(&error, prec_internal);
    arpra_init2(&zNew, z->precision);
    mpfr_set_ui(&error, 0, MPFR_RNDU);
    mpfr_set_ui(&(zNew.radius), 0, MPFR_RNDU);

    // z_0 = (alpha * x_0) + gamma
    if (mpfr_fma(&(zNew.centre), alpha, &(x->centre), gamma, MPFR_RNDN)) {
        arpra_helper_error_half_ulp(&temp1, &(zNew.centre));
        mpfr_add(&error, &error, &temp1, MPFR_RNDU);
    }

    // Allocate memory for all possible deviation terms.
    zNew.nTerms = x->nTerms + 1;
    zNew.symbols = malloc(zNew.nTerms * sizeof(arpra_uint));
    zNew.deviations = malloc(zNew.nTerms * sizeof(arpra_mpfr));

    for (xTerm = 0, zTerm = 0; xTerm < x->nTerms; xTerm++) {
        zNew.symbols[zTerm] = x->symbols[xTerm];
        mpfr_init2(&(zNew.deviations[zTerm]), prec_internal);

        // z_i = (alpha * x_i)
        if (mpfr_mul(&(zNew.deviations[zTerm]), alpha, &(x->deviations[xTerm]), MPFR_RNDN)) {
            arpra_helper_error_half_ulp(&temp1, &(zNew.deviations[zTerm]));
            mpfr_add(&error, &error, &temp1, MPFR_RNDU);
        }

        mpfr_abs(&temp1, &(zNew.deviations[zTerm]), MPFR_RNDU);
        mpfr_add(&(zNew.radius), &(zNew.radius), &temp1, MPFR_RNDU);
        zTerm++;
    }

    // Add delta error.
    mpfr_add(&error, &error, delta, MPFR_RNDU);

    // Store numerical error term.
    zNew.symbols[zTerm] = arpra_next_symbol();
    zNew.deviations[zTerm] = error;
    mpfr_add(&(zNew.radius), &(zNew.radius), &(zNew.deviations[zTerm]), MPFR_RNDU);
    zNew.nTerms = zTerm + 1;

    // Compute true_range in working precision.
    arpra_helper_true_range(&zNew);

    // Handle domain violations.
    if (mpfr_nan_p(&(zNew.centre)) || mpfr_nan_p(&(zNew.radius))) {
        arpra_set_nan(&zNew);
    }
    else if (mpfr_inf_p(&(zNew.centre)) || mpfr_inf_p(&(zNew.radius))) {
        arpra_set_inf(&zNew);
    }

    // Clear vars, and set z.
    mpfr_clear(&temp1);
    mpfr_clear(&temp2);
    arpra_clear(z);
    *z = zNew;
}
