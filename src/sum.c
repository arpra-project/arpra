/*
 * sum.c -- Sum an array of Arpra ranges.
 *
 * Copyright 2016-2020 James Paul Turner.
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

void arpra_sum_exact (arpra_range *y, arpra_range *x, const arpra_uint n)
{
    arpra_uint i, j;
    arpra_uint xSymbol, yTerm;
    arpra_uint *xTerm;
    arpra_int xHasNext;
    arpra_mpfr **summands;
    arpra_mpfr temp1, temp2, error;
    arpra_range yy;
    arpra_prec prec_internal;

    // Domain violations:
    // (NaN) + ... + (NaN) = (NaN)
    // (NaN) + ... + (R)   = (NaN)
    // (Inf) + ... + (Inf) = (NaN)
    // (Inf) + ... + (R)   = (Inf)

    // Handle domain violations.
    for (i = 0; i < n; i++) {
        if (arpra_nan_p(&x[i])) {
            arpra_set_nan(y);
            return;
        }
    }
    for (i = 0; i < n; i++) {
        if (arpra_inf_p(&x[i])) {
            for (++i; i < n; i++) {
                if (arpra_inf_p(&x[i])) {
                    arpra_set_nan(y);
                    return;
                }
            }
            arpra_set_inf(y);
            return;
        }
    }

    // Handle trivial cases.
    if (n <= 2) {
        if (n == 2) {
            arpra_add(y, &x[0], &x[1]);
            return;
        }
        else if (n == 1) {
            arpra_set(y, &x[0]);
            return;
        }
        else {
            arpra_set_zero(y);
            return;
        }
    }

    // Initialise vars.
    prec_internal = arpra_get_internal_precision();
    mpfr_init2(&temp1, prec_internal + 8);
    mpfr_init2(&temp2, prec_internal + 8);
    mpfr_init2(&error, prec_internal);
    arpra_init2(&yy, y->precision);
    mpfr_set_zero(&error, 1);
    mpfr_set_zero(&(yy.radius), 1);
    xTerm = malloc(n * sizeof(arpra_uint));
    summands = malloc(n * sizeof(arpra_mpfr *));

    // Zero term indexes, and fill summand array with centre values.
    yTerm = 0;
    for (i = 0; i < n; i++) {
        xTerm[i] = 0;
        summands[i] = &(x[i].centre);
    }

    // y[0] = x1[0] + ... + xn[0]
    if (mpfr_sum(&(yy.centre), summands, n, MPFR_RNDN)) {
        arpra_helper_error_half_ulp(&temp1, &(yy.centre));
        mpfr_add(&error, &error, &temp1, MPFR_RNDU);
    }

    // Allocate memory for all possible deviation terms.
    yy.nTerms = 1;
    for (i = 0; i < n; i++) {
        yy.nTerms += x[i].nTerms;
    }
    yy.symbols = malloc(yy.nTerms * sizeof(arpra_uint));
    yy.deviations = malloc(yy.nTerms * sizeof(arpra_mpfr));

    // For all unique symbols in x.
    xHasNext = yy.nTerms > 1;
    while (xHasNext) {
        xHasNext = 0;
        xSymbol = -1;

        // Find and set the next symbol in y.
        for (i = 0; i < n; i++) {
            if (xTerm[i] < x[i].nTerms) {
                if (x[i].symbols[xTerm[i]] < xSymbol) {
                    xSymbol = x[i].symbols[xTerm[i]];
                }
            }
        }
        yy.symbols[yTerm] = xSymbol;
        mpfr_init2(&(yy.deviations[yTerm]), prec_internal);

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

        // y[i] = x1[i] + ... + xn[i]
        if (mpfr_sum(&(yy.deviations[yTerm]), summands, j, MPFR_RNDN)) {
            arpra_helper_error_half_ulp(&temp1, &(yy.deviations[yTerm]));
            mpfr_add(&error, &error, &temp1, MPFR_RNDU);
        }

        mpfr_abs(&temp1, &(yy.deviations[yTerm]), MPFR_RNDU);
        mpfr_add(&(yy.radius), &(yy.radius), &temp1, MPFR_RNDU);
        yTerm++;
    }

    // Store numerical error term.
    yy.symbols[yTerm] = arpra_helper_next_symbol();
    yy.deviations[yTerm] = error;
    mpfr_add(&(yy.radius), &(yy.radius), &(yy.deviations[yTerm]), MPFR_RNDU);
    yy.nTerms = yTerm + 1;

    // Compute true_range, and add rounding error.
    arpra_helper_range_rounded(&yy);

    // Check for NaN and Inf.
    arpra_helper_check_result(&yy);

    // Clear vars, and set y.
    mpfr_clear(&temp1);
    mpfr_clear(&temp2);
    arpra_clear(y);
    *y = yy;
    free(xTerm);
    free(summands);
}

void arpra_sum_recursive (arpra_range *y, arpra_range *x, const arpra_uint n)
{
    arpra_uint i, j;
    arpra_uint xSymbol, yTerm;
    arpra_uint *xTerm;
    arpra_int xHasNext;
    arpra_mpfr **summands;
    arpra_mpfr temp1, temp2, error;
    arpra_range yy;
    arpra_prec prec_internal;

    // Domain violations:
    // (NaN) + ... + (NaN) = (NaN)
    // (NaN) + ... + (R)   = (NaN)
    // (Inf) + ... + (Inf) = (NaN)
    // (Inf) + ... + (R)   = (Inf)

    // Handle domain violations.
    for (i = 0; i < n; i++) {
        if (arpra_nan_p(&x[i])) {
            arpra_set_nan(y);
            return;
        }
    }
    for (i = 0; i < n; i++) {
        if (arpra_inf_p(&x[i])) {
            for (++i; i < n; i++) {
                if (arpra_inf_p(&x[i])) {
                    arpra_set_nan(y);
                    return;
                }
            }
            arpra_set_inf(y);
            return;
        }
    }

    // Handle trivial cases.
    if (n <= 2) {
        if (n == 2) {
            arpra_add(y, &x[0], &x[1]);
            return;
        }
        else if (n == 1) {
            arpra_set(y, &x[0]);
            return;
        }
        else {
            arpra_set_zero(y);
            return;
        }
    }

    // Initialise vars.
    prec_internal = arpra_get_internal_precision();
    mpfr_init2(&temp1, prec_internal + 8);
    mpfr_init2(&temp2, prec_internal + 8);
    mpfr_init2(&error, prec_internal);
    arpra_init2(&yy, y->precision);
    mpfr_set_zero(&error, 1);
    mpfr_set_zero(&(yy.radius), 1);
    xTerm = malloc(n * sizeof(arpra_uint));
    summands = malloc(n * sizeof(arpra_mpfr *));

    // Zero term indexes, and fill summand array with centre values.
    yTerm = 0;
    for (i = 0; i < n; i++) {
        xTerm[i] = 0;
        summands[i] = &(x[i].centre);
    }

    // y[0] = x1[0] + ... + xn[0]
    if (mpfr_sum(&(yy.centre), summands, n, MPFR_RNDN)) {
        arpra_helper_error_half_ulp(&temp1, &(yy.centre));
        mpfr_add(&error, &error, &temp1, MPFR_RNDU);
    }

    // Allocate memory for all possible deviation terms.
    yy.nTerms = 1;
    for (i = 0; i < n; i++) {
        yy.nTerms += x[i].nTerms;
    }
    yy.symbols = malloc(yy.nTerms * sizeof(arpra_uint));
    yy.deviations = malloc(yy.nTerms * sizeof(arpra_mpfr));

    // For all unique symbols in x.
    xHasNext = yy.nTerms > 1;
    while (xHasNext) {
        xHasNext = 0;
        xSymbol = -1;

        // Find and set the next symbol in y.
        for (i = 0; i < n; i++) {
            if (xTerm[i] < x[i].nTerms) {
                if (x[i].symbols[xTerm[i]] < xSymbol) {
                    xSymbol = x[i].symbols[xTerm[i]];
                }
            }
        }
        yy.symbols[yTerm] = xSymbol;
        mpfr_init2(&(yy.deviations[yTerm]), prec_internal);

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

        // y[i] = x1[i] + ... + xn[i]
        if (mpfr_sum(&(yy.deviations[yTerm]), summands, j, MPFR_RNDN)) {
            arpra_helper_error_half_ulp(&temp1, &(yy.deviations[yTerm]));
            mpfr_add(&error, &error, &temp1, MPFR_RNDU);
        }

        mpfr_abs(&temp1, &(yy.deviations[yTerm]), MPFR_RNDU);
        mpfr_add(&(yy.radius), &(yy.radius), &temp1, MPFR_RNDU);
        yTerm++;
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

        if (mpfr_sgn(&(x[i].centre)) >= 0) {
            mpfr_add(&(sum_error[i]), &(x[i].centre), &(x[i].radius), MPFR_RNDU);
        }
        else {
            mpfr_sub(&(sum_error[i]), &(x[i].centre), &(x[i].radius), MPFR_RNDD);
            mpfr_abs(&(sum_error[i]), &(sum_error[i]), MPFR_RNDU);
        }
    }

    // Compute error(sum(x)) = (n - 1) u sum(|x|).
    mpfr_set_si_2exp(&temp1, 1, -yy.precision, MPFR_RNDU);
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

    // Store numerical error term.
    yy.symbols[yTerm] = arpra_helper_next_symbol();
    yy.deviations[yTerm] = error;
    mpfr_add(&(yy.radius), &(yy.radius), &(yy.deviations[yTerm]), MPFR_RNDU);
    yy.nTerms = yTerm + 1;

    // Compute true_range, and add rounding error.
    arpra_helper_range_rounded(&yy);

    // Check for NaN and Inf.
    arpra_helper_check_result(&yy);

    // Clear vars, and set y.
    mpfr_clear(&temp1);
    mpfr_clear(&temp2);
    arpra_clear(y);
    *y = yy;
    free(xTerm);
    free(summands);
}
