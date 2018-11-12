/*
 * affine_2.c -- Compute a bivariate affine function of two Arpra ranges.
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

void arpra_affine_2 (arpra_range *z, const arpra_range *x, const arpra_range *y,
                     const arpra_mpfr *alpha, const arpra_mpfr *beta,
                     const arpra_mpfr *gamma, const arpra_mpfr *delta)
{
    arpra_uint xTerm, yTerm, zTerm;
    arpra_int xHasNext, yHasNext;
    arpra_mpfr z_lo, z_hi, temp, error;
    arpra_prec prec_internal;
    arpra_range zNew;

    // Domain violations:
    // NaN  +  NaN  =  NaN
    // NaN  +  r    =  NaN
    // Inf  +  Inf  =  NaN
    // Inf  +  r    =  Inf
    // s.t. (r in R)

    // Handle domain violations.
    if (arpra_nan_p(x) || arpra_nan_p(y)) {
        arpra_set_nan(z);
        return;
    }
    if (arpra_inf_p(x)) {
        if (arpra_inf_p(y)) {
            arpra_set_nan(z);
        }
        else {
            arpra_set_inf(z);
        }
        return;
    }
    if (arpra_inf_p(y)) {
        if (arpra_inf_p(x)) {
            arpra_set_nan(z);
        }
        else {
            arpra_set_inf(z);
        }
        return;
    }

    // Initialise vars.
    prec_internal = arpra_get_internal_precision();
    mpfr_init2(&z_lo, prec_internal);
    mpfr_init2(&z_hi, prec_internal);
    mpfr_init2(&temp, prec_internal);
    mpfr_init2(&error, prec_internal);
    mpfr_init2(&(zNew.centre), prec_internal);
    mpfr_init2(&(zNew.radius), prec_internal);
    mpfr_set(&error, delta, MPFR_RNDU);
    mpfr_set_si(&(zNew.radius), 0, MPFR_RNDU);
    zNew.precision = z->precision;

    // z_0 = (alpha * x_0) + (beta * y_0) + gamma
    if (arpra_helper_term(&(zNew.centre), &(x->centre), &(y->centre), alpha, beta, gamma)) {
        arpra_helper_error_half_ulp(&temp, &(zNew.centre));
        mpfr_add(&error, &error, &temp, MPFR_RNDU);
    }

    // Allocate memory for all possible deviation terms.
    zNew.nTerms = x->nTerms + y->nTerms + 1;
    zNew.symbols = malloc(zNew.nTerms * sizeof(arpra_uint));
    zNew.deviations = malloc(zNew.nTerms * sizeof(arpra_mpfr));

    xTerm = 0;
    yTerm = 0;
    zTerm = 0;
    xHasNext = x->nTerms > 0;
    yHasNext = y->nTerms > 0;
    while (xHasNext || yHasNext) {
        if ((!yHasNext) || (xHasNext && (x->symbols[xTerm] < y->symbols[yTerm]))) {
            zNew.symbols[zTerm] = x->symbols[xTerm];
            mpfr_init2(&(zNew.deviations[zTerm]), prec_internal);

            // z_i = (alpha * x_i)
            if (mpfr_mul(&(zNew.deviations[zTerm]), alpha, &(x->deviations[xTerm]), MPFR_RNDN)) {
                arpra_helper_error_half_ulp(&temp, &(zNew.deviations[zTerm]));
                mpfr_add(&error, &error, &temp, MPFR_RNDU);
            }

            xHasNext = ++xTerm < x->nTerms;
        }
        else if ((!xHasNext) || (yHasNext && (y->symbols[yTerm] < x->symbols[xTerm]))) {
            zNew.symbols[zTerm] = y->symbols[yTerm];
            mpfr_init2(&(zNew.deviations[zTerm]), prec_internal);

            // z_i = (beta * y_i)
            if (mpfr_mul(&(zNew.deviations[zTerm]), beta, &(y->deviations[yTerm]), MPFR_RNDN)) {
                arpra_helper_error_half_ulp(&temp, &(zNew.deviations[zTerm]));
                mpfr_add(&error, &error, &temp, MPFR_RNDU);
            }

            yHasNext = ++yTerm < y->nTerms;
        }
        else {
            zNew.symbols[zTerm] = x->symbols[xTerm];
            mpfr_init2(&(zNew.deviations[zTerm]), prec_internal);

            // z_i = (alpha * x_i) + (beta * y_i)
            if (arpra_helper_term(&(zNew.deviations[zTerm]), &(x->deviations[xTerm]), &(y->deviations[yTerm]), alpha, beta, NULL)) {
                arpra_helper_error_half_ulp(&temp, &(zNew.deviations[zTerm]));
                mpfr_add(&error, &error, &temp, MPFR_RNDU);
            }

            xHasNext = ++xTerm < x->nTerms;
            yHasNext = ++yTerm < y->nTerms;
        }

        // Store nonzero deviation terms.
        if (mpfr_zero_p(&(zNew.deviations[zTerm]))) {
            mpfr_clear(&(zNew.deviations[zTerm]));
        }
        else {
            mpfr_abs(&temp, &(zNew.deviations[zTerm]), MPFR_RNDU);
            mpfr_add(&(zNew.radius), &(zNew.radius), &temp, MPFR_RNDU);
            zTerm++;
        }
    }

    // Round the result to the target precision.
    mpfr_sub(&z_lo, &(zNew.centre), &(zNew.radius), MPFR_RNDD);
    mpfr_sub(&z_lo, &z_lo, &error, MPFR_RNDD);
    if (mpfr_prec_round(&z_lo, zNew.precision, MPFR_RNDD)) {
        arpra_helper_error_ulp(&z_lo, &z_lo);
    }
    else {
        mpfr_set_ui(&z_lo, 0, MPFR_RNDN);
    }
    mpfr_add(&z_hi, &(zNew.centre), &(zNew.radius), MPFR_RNDU);
    mpfr_add(&z_hi, &z_hi, &error, MPFR_RNDU);
    if (mpfr_prec_round(&z_hi, zNew.precision, MPFR_RNDU)) {
        arpra_helper_error_ulp(&z_hi, &z_hi);
    }
    else {
        mpfr_set_ui(&z_hi, 0, MPFR_RNDN);
    }
    mpfr_max(&temp, &z_lo, &z_hi, MPFR_RNDU);
    mpfr_add(&error, &error, &temp, MPFR_RNDU);

    // Store nonzero numerical error term.
    if (!mpfr_zero_p(&error)) {
        zNew.symbols[zTerm] = arpra_next_symbol();
        mpfr_init2(&(zNew.deviations[zTerm]), prec_internal);
        mpfr_set(&(zNew.deviations[zTerm]), &error, MPFR_RNDU);
        mpfr_add(&(zNew.radius), &(zNew.radius), &(zNew.deviations[zTerm]), MPFR_RNDU);
        zTerm++;
    }

    // Handle domain violations, and free unused memory.
    zNew.nTerms = zTerm;
    if (mpfr_nan_p(&(zNew.centre)) || mpfr_nan_p(&(zNew.radius))) {
        arpra_set_nan(&zNew);
    }
    else if (mpfr_inf_p(&(zNew.centre)) || mpfr_inf_p(&(zNew.radius))) {
        arpra_set_inf(&zNew);
    }
    else {
        if (zNew.nTerms == 0) {
            free(zNew.symbols);
            free(zNew.deviations);
        }
    }

    // Clear vars, and set z.
    mpfr_clear(&z_lo);
    mpfr_clear(&z_hi);
    mpfr_clear(&temp);
    mpfr_clear(&error);
    arpra_clear(z);
    *z = zNew;
}
