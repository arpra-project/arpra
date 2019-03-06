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
    arpra_mpfr temp1, temp2, error;
    arpra_range zNew;
    arpra_prec prec_internal;

    // Domain violations:
    // NaN  +  NaN  +  ...  =  NaN
    // NaN  +  r    +  ...  =  NaN
    // Inf  +  Inf  +  ...  =  NaN
    // Inf  +  r    +  ...  =  Inf

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
    mpfr_init2(&temp1, prec_internal + 8);
    mpfr_init2(&temp2, prec_internal + 8);
    mpfr_init2(&error, prec_internal);
    arpra_init2(&zNew, z->precision);
    mpfr_set_ui(&error, 0, MPFR_RNDU);
    mpfr_set_ui(&(zNew.radius), 0, MPFR_RNDU);
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
        arpra_helper_error_half_ulp(&temp1, &(zNew.centre));
        mpfr_add(&error, &error, &temp1, MPFR_RNDU);
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
            arpra_helper_error_half_ulp(&temp1, &(zNew.deviations[zTerm]));
            mpfr_add(&error, &error, &temp1, MPFR_RNDU);
        }

        mpfr_abs(&temp1, &(zNew.deviations[zTerm]), MPFR_RNDU);
        mpfr_add(&(zNew.radius), &(zNew.radius), &temp1, MPFR_RNDU);
        zTerm++;
    }

    // Round range to target precision.
    mpfr_sub(&temp1, &(zNew.centre), &(zNew.radius), MPFR_RNDD);
    mpfr_sub(&temp1, &temp1, &error, MPFR_RNDD);
    mpfr_set(&(zNew.true_range.left), &temp1, MPFR_RNDD);
    mpfr_sub(&temp1, &temp1, &(zNew.true_range.left), MPFR_RNDU);
    mpfr_add(&temp2, &(zNew.centre), &(zNew.radius), MPFR_RNDU);
    mpfr_add(&temp2, &temp2, &error, MPFR_RNDU);
    mpfr_set(&(zNew.true_range.right), &temp2, MPFR_RNDU);
    mpfr_sub(&temp2, &(zNew.true_range.right), &temp2, MPFR_RNDU);
    mpfr_max(&temp1, &temp1, &temp2, MPFR_RNDU);
    mpfr_add(&error, &error, &temp1, MPFR_RNDU);

    // Store numerical error term.
    zNew.symbols[zTerm] = arpra_next_symbol();
    zNew.deviations[zTerm] = error;
    mpfr_add(&(zNew.radius), &(zNew.radius), &(zNew.deviations[zTerm]), MPFR_RNDU);
    zNew.nTerms = zTerm + 1;

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
    arpra_mpfr temp1, temp2, error;
    arpra_range zNew;
    arpra_prec prec_internal;

    // Domain violations:
    // NaN  +  NaN  +  ...  =  NaN
    // NaN  +  r    +  ...  =  NaN
    // Inf  +  Inf  +  ...  =  NaN
    // Inf  +  r    +  ...  =  Inf

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
    mpfr_init2(&temp1, prec_internal + 8);
    mpfr_init2(&temp2, prec_internal + 8);
    mpfr_init2(&error, prec_internal);
    arpra_init2(&zNew, z->precision);
    mpfr_set_ui(&error, 0, MPFR_RNDU);
    mpfr_set_ui(&(zNew.radius), 0, MPFR_RNDU);
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
        arpra_helper_error_half_ulp(&temp1, &(zNew.centre));
        mpfr_add(&error, &error, &temp1, MPFR_RNDU);
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
            arpra_helper_error_half_ulp(&temp1, &(zNew.deviations[zTerm]));
            mpfr_add(&error, &error, &temp1, MPFR_RNDU);
        }

        mpfr_abs(&temp1, &(zNew.deviations[zTerm]), MPFR_RNDU);
        mpfr_add(&(zNew.radius), &(zNew.radius), &temp1, MPFR_RNDU);
        zTerm++;
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
    mpfr_set_si_2exp(&temp1, 1, -zNew.precision, MPFR_RNDU);
    mpfr_mul_ui(&temp1, &temp1, (n - 1), MPFR_RNDU);
    arpra_helper_mpfr_sum(&(sum_error[0]), sum_error, n, MPFR_RNDU);
    mpfr_mul(&(sum_error[0]), &(sum_error[0]), &temp1, MPFR_RNDU);
    mpfr_add(&error, &error, &(sum_error[0]), MPFR_RNDU);

    for (i = 0; i < n; i++) {
        mpfr_clear(&(sum_error[i]));
    }

    /*
     * END: Error bound for recursive summation.
     */

    // Round range to target precision.
    mpfr_sub(&temp1, &(zNew.centre), &(zNew.radius), MPFR_RNDD);
    mpfr_sub(&temp1, &temp1, &error, MPFR_RNDD);
    mpfr_set(&(zNew.true_range.left), &temp1, MPFR_RNDD);
    mpfr_sub(&temp1, &temp1, &(zNew.true_range.left), MPFR_RNDU);
    mpfr_add(&temp2, &(zNew.centre), &(zNew.radius), MPFR_RNDU);
    mpfr_add(&temp2, &temp2, &error, MPFR_RNDU);
    mpfr_set(&(zNew.true_range.right), &temp2, MPFR_RNDU);
    mpfr_sub(&temp2, &(zNew.true_range.right), &temp2, MPFR_RNDU);
    mpfr_max(&temp1, &temp1, &temp2, MPFR_RNDU);
    mpfr_add(&error, &error, &temp1, MPFR_RNDU);

    // Store numerical error term.
    zNew.symbols[zTerm] = arpra_next_symbol();
    zNew.deviations[zTerm] = error;
    mpfr_add(&(zNew.radius), &(zNew.radius), &(zNew.deviations[zTerm]), MPFR_RNDU);
    zNew.nTerms = zTerm + 1;

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
    free(xTerm);
    free(summands);
}
