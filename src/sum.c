/*
 * sum.c -- Sum an array of Arpra ranges.
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

void arpra_sum_exact (arpra_range *z, arpra_range *x, const arpra_uint n)
{
    arpra_uint i, j;
    arpra_uint xSymbol, zTerm;
    arpra_uint *xTerm;
    arpra_int xHasNext;
    arpra_mpfr **summands;
    arpra_mpfr temp, error;
    arpra_mpfi z_range;
    arpra_range zNew;
    arpra_prec prec_internal;

    // Domain violations:
    // NaN  +  NaN  +  ...  =  NaN
    // NaN  +  r    +  ...  =  NaN
    // Inf  +  Inf  +  ...  =  NaN
    // Inf  +  r    +  ...  =  Inf
    // s.t. (r in R)

    // Handle domain violations.
    for (i = 0; i < n; i++) {
        if (arpra_nan_p(&x[i])) {
            arpra_set_nan(z);
            return;
        }
    }
    for (i = 0; i < n; i++) {
        if (arpra_inf_p(&x[i])) {
            for (++i; i < n; i++) {
                if (arpra_inf_p(&x[i])) {
                    arpra_set_nan(z);
                    return;
                }
            }
            arpra_set_inf(z);
            return;
        }
    }

    // Handle trivial cases.
    if (n <= 2) {
        if (n == 2) {
            arpra_add(z, &x[0], &x[1]);
            return;
        }
        else if (n == 1) {
            arpra_set(z, &x[0]);
            return;
        }
        else {
            arpra_set_zero(z);
            return;
        }
    }

    // Initialise vars.
    prec_internal = arpra_get_internal_precision();
    mpfr_init2(&temp, prec_internal);
    mpfr_init2(&error, prec_internal);
    mpfi_init2(&z_range, prec_internal);
    arpra_init2(&zNew, z->precision);
    mpfr_set_si(&error, 0, MPFR_RNDU);
    mpfr_set_si(&(zNew.radius), 0, MPFR_RNDU);
    xTerm = malloc(n * sizeof(arpra_uint));
    summands = malloc(n * sizeof(arpra_mpfr *));

    // Zero term indexes, and fill summand array with centre values.
    zTerm = 0;
    for (i = 0; i < n; i++) {
        xTerm[i] = 0;
        summands[i] = &(x[i].centre);
    }

    // z_0 = x[1]_0 + ... + x[n]_0
    if (mpfr_sum(&(zNew.centre), summands, n, MPFR_RNDN)) {
        arpra_helper_error_half_ulp(&temp, &(zNew.centre));
        mpfr_add(&error, &error, &temp, MPFR_RNDU);
    }

    // Allocate memory for all possible deviation terms.
    zNew.nTerms = 1;
    for (i = 0; i < n; i++) {
        zNew.nTerms += x[i].nTerms;
    }
    zNew.symbols = malloc(zNew.nTerms * sizeof(arpra_uint));
    zNew.deviations = malloc(zNew.nTerms * sizeof(arpra_mpfr));

    // For all unique symbols in x.
    xHasNext = zNew.nTerms > 1;
    while (xHasNext) {
        xHasNext = 0;
        xSymbol = -1;

        // Find and set the next symbol in z.
        for (i = 0; i < n; i++) {
            if (xTerm[i] < x[i].nTerms) {
                if (x[i].symbols[xTerm[i]] < xSymbol) {
                    xSymbol = x[i].symbols[xTerm[i]];
                }
            }
        }
        zNew.symbols[zTerm] = xSymbol;
        mpfr_init2(&(zNew.deviations[zTerm]), prec_internal);

        // For all x with the next symbol:
        for (i = 0, j = 0; i < n; i++) {
            if (x[i].symbols[xTerm[i]] == xSymbol) {
                // Get next deviation pointer of x[i].
                summands[j++] = &(x[i].deviations[xTerm[i]]);

                // Check for more symbols in x[i].
                if (++xTerm[i] < x[i].nTerms) {
                    xHasNext = 1;
                }
            }
        }

        // z_i = x[1]_i + ... + x[n]_i
        if (mpfr_sum(&(zNew.deviations[zTerm]), summands, j, MPFR_RNDN)) {
            arpra_helper_error_half_ulp(&temp, &(zNew.deviations[zTerm]));
            mpfr_add(&error, &error, &temp, MPFR_RNDU);
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
    mpfr_sub(&(z_range.left), &(zNew.centre), &(zNew.radius), MPFR_RNDD);
    mpfr_sub(&(z_range.left), &(z_range.left), &error, MPFR_RNDD);
    if (mpfr_prec_round(&(z_range.left), zNew.precision, MPFR_RNDD)) {
        arpra_helper_error_ulp(&(z_range.left), &(z_range.left));
    }
    else {
        mpfr_set_ui(&(z_range.left), 0, MPFR_RNDN);
    }
    mpfr_add(&(z_range.right), &(zNew.centre), &(zNew.radius), MPFR_RNDU);
    mpfr_add(&(z_range.right), &(z_range.right), &error, MPFR_RNDU);
    if (mpfr_prec_round(&(z_range.right), zNew.precision, MPFR_RNDU)) {
        arpra_helper_error_ulp(&(z_range.right), &(z_range.right));
    }
    else {
        mpfr_set_ui(&(z_range.right), 0, MPFR_RNDN);
    }
    mpfr_max(&temp, &(z_range.left), &(z_range.right), MPFR_RNDU);
    mpfr_add(&error, &error, &temp, MPFR_RNDU);

    // Store nonzero numerical error term.
    if (!mpfr_zero_p(&error)) {
        zNew.symbols[zTerm] = arpra_next_symbol();
        mpfr_init2(&(zNew.deviations[zTerm]), prec_internal);
        mpfr_set(&(zNew.deviations[zTerm]), &error, MPFR_RNDU);
        mpfr_add(&(zNew.radius), &(zNew.radius), &(zNew.deviations[zTerm]), MPFR_RNDU);
        zTerm++;
    }

    // Compute true range.
    mpfr_sub(&(zNew.true_range.left), &(zNew.centre), &(zNew.radius), MPFR_RNDD);
    mpfr_add(&(zNew.true_range.right), &(zNew.centre), &(zNew.radius), MPFR_RNDU);

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
    mpfr_clear(&temp);
    mpfr_clear(&error);
    mpfi_clear(&z_range);
    arpra_clear(z);
    *z = zNew;
    free(xTerm);
    free(summands);
}

void arpra_sum_recursive (arpra_range *z, arpra_range *x, const arpra_uint n)
{
    arpra_uint i, j;
    arpra_uint xSymbol, zTerm;
    arpra_uint *xTerm;
    arpra_int xHasNext;
    arpra_mpfr **summands;
    arpra_mpfr temp, error;
    arpra_mpfi z_range;
    arpra_range zNew;
    arpra_prec prec_internal;

    // Domain violations:
    // NaN  +  NaN  +  ...  =  NaN
    // NaN  +  r    +  ...  =  NaN
    // Inf  +  Inf  +  ...  =  NaN
    // Inf  +  r    +  ...  =  Inf
    // s.t. (r in R)

    // Handle domain violations.
    for (i = 0; i < n; i++) {
        if (arpra_nan_p(&x[i])) {
            arpra_set_nan(z);
            return;
        }
    }
    for (i = 0; i < n; i++) {
        if (arpra_inf_p(&x[i])) {
            for (++i; i < n; i++) {
                if (arpra_inf_p(&x[i])) {
                    arpra_set_nan(z);
                    return;
                }
            }
            arpra_set_inf(z);
            return;
        }
    }

    // Handle trivial cases.
    if (n <= 2) {
        if (n == 2) {
            arpra_add(z, &x[0], &x[1]);
            return;
        }
        else if (n == 1) {
            arpra_set(z, &x[0]);
            return;
        }
        else {
            arpra_set_zero(z);
            return;
        }
    }

    // Initialise vars.
    prec_internal = arpra_get_internal_precision();
    mpfr_init2(&temp, prec_internal);
    mpfr_init2(&error, prec_internal);
    mpfi_init2(&z_range, prec_internal);
    arpra_init2(&zNew, z->precision);
    mpfr_set_si(&error, 0, MPFR_RNDU);
    mpfr_set_si(&(zNew.radius), 0, MPFR_RNDU);
    xTerm = malloc(n * sizeof(arpra_uint));
    summands = malloc(n * sizeof(arpra_mpfr *));

    // Zero term indexes, and fill summand array with centre values.
    zTerm = 0;
    for (i = 0; i < n; i++) {
        xTerm[i] = 0;
        summands[i] = &(x[i].centre);
    }

    // z_0 = x[1]_0 + ... + x[n]_0
    if (mpfr_sum(&(zNew.centre), summands, n, MPFR_RNDN)) {
        arpra_helper_error_half_ulp(&temp, &(zNew.centre));
        mpfr_add(&error, &error, &temp, MPFR_RNDU);
    }

    // Allocate memory for all possible deviation terms.
    zNew.nTerms = 1;
    for (i = 0; i < n; i++) {
        zNew.nTerms += x[i].nTerms;
    }
    zNew.symbols = malloc(zNew.nTerms * sizeof(arpra_uint));
    zNew.deviations = malloc(zNew.nTerms * sizeof(arpra_mpfr));

    // For all unique symbols in x.
    xHasNext = zNew.nTerms > 1;
    while (xHasNext) {
        xHasNext = 0;
        xSymbol = -1;

        // Find and set the next symbol in z.
        for (i = 0; i < n; i++) {
            if (xTerm[i] < x[i].nTerms) {
                if (x[i].symbols[xTerm[i]] < xSymbol) {
                    xSymbol = x[i].symbols[xTerm[i]];
                }
            }
        }
        zNew.symbols[zTerm] = xSymbol;
        mpfr_init2(&(zNew.deviations[zTerm]), prec_internal);

        // For all x with the next symbol:
        for (i = 0, j = 0; i < n; i++) {
            if (x[i].symbols[xTerm[i]] == xSymbol) {
                // Get next deviation pointer of x[i].
                summands[j++] = &(x[i].deviations[xTerm[i]]);

                // Check for more symbols in x[i].
                if (++xTerm[i] < x[i].nTerms) {
                    xHasNext = 1;
                }
            }
        }

        // z_i = x[1]_i + ... + x[n]_i
        if (mpfr_sum(&(zNew.deviations[zTerm]), summands, j, MPFR_RNDN)) {
            arpra_helper_error_half_ulp(&temp, &(zNew.deviations[zTerm]));
            mpfr_add(&error, &error, &temp, MPFR_RNDU);
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

    /*
     * Error bound for recursive summation (any ordering, any n).
     *
     * Rump, S. Error estimation of floating-point summation and dot product.
     * BIT Numer Math (2012) 52:201-220.
     */

    // Compute |x|.
    arpra_mpfr *sum_error;
    sum_error = arpra_helper_buffer_mpfr(n);
    for (i = 0; i < n; i++) {
        mpfr_init2(&(sum_error[i]), prec_internal);

        if (mpfr_sgn(&(x->centre)) >= 0) {
            mpfr_add(&(sum_error[i]), &(x->centre), &(x->radius), MPFR_RNDU);
        }
        else {
            mpfr_sub(&(sum_error[i]), &(x->centre), &(x->radius), MPFR_RNDD);
            mpfr_abs(&(sum_error[i]), &(sum_error[i]), MPFR_RNDU);
        }
    }

    // Compute error(sum(x)) = (n - 1) u sum(|x|).
    mpfr_set_si_2exp(&temp, 1, -zNew.precision, MPFR_RNDU);
    mpfr_mul_ui(&temp, &temp, (n - 1), MPFR_RNDU);
    arpra_helper_mpfr_sum(&(sum_error[0]), sum_error, n, MPFR_RNDU);
    mpfr_mul(&(sum_error[0]), &(sum_error[0]), &temp, MPFR_RNDU);
    mpfr_add(&error, &error, &(sum_error[0]), MPFR_RNDU);

    for (i = 0; i < n; i++) {
        mpfr_clear(&(sum_error[i]));
    }

    /*
     * END: Error bound for recursive summation.
     */

    // Round the result to the target precision.
    mpfr_sub(&(z_range.left), &(zNew.centre), &(zNew.radius), MPFR_RNDD);
    mpfr_sub(&(z_range.left), &(z_range.left), &error, MPFR_RNDD);
    if (mpfr_prec_round(&(z_range.left), zNew.precision, MPFR_RNDD)) {
        arpra_helper_error_ulp(&(z_range.left), &(z_range.left));
    }
    else {
        mpfr_set_ui(&(z_range.left), 0, MPFR_RNDN);
    }
    mpfr_add(&(z_range.right), &(zNew.centre), &(zNew.radius), MPFR_RNDU);
    mpfr_add(&(z_range.right), &(z_range.right), &error, MPFR_RNDU);
    if (mpfr_prec_round(&(z_range.right), zNew.precision, MPFR_RNDU)) {
        arpra_helper_error_ulp(&(z_range.right), &(z_range.right));
    }
    else {
        mpfr_set_ui(&(z_range.right), 0, MPFR_RNDN);
    }
    mpfr_max(&temp, &(z_range.left), &(z_range.right), MPFR_RNDU);
    mpfr_add(&error, &error, &temp, MPFR_RNDU);

    // Store nonzero numerical error term.
    if (!mpfr_zero_p(&error)) {
        zNew.symbols[zTerm] = arpra_next_symbol();
        mpfr_init2(&(zNew.deviations[zTerm]), prec_internal);
        mpfr_set(&(zNew.deviations[zTerm]), &error, MPFR_RNDU);
        mpfr_add(&(zNew.radius), &(zNew.radius), &(zNew.deviations[zTerm]), MPFR_RNDU);
        zTerm++;
    }

    // Compute true range.
    mpfr_sub(&(zNew.true_range.left), &(zNew.centre), &(zNew.radius), MPFR_RNDD);
    mpfr_add(&(zNew.true_range.right), &(zNew.centre), &(zNew.radius), MPFR_RNDU);

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
    mpfr_clear(&temp);
    mpfr_clear(&error);
    mpfi_clear(&z_range);
    arpra_clear(z);
    *z = zNew;
    free(xTerm);
    free(summands);
}
